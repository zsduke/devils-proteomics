"""Local table analysis and coordinate-viewer entry point."""
from __future__ import annotations
import argparse
import base64
import binascii
import hmac
import importlib.util
import io
import json
import mimetypes
from pathlib import Path
import re
import secrets
import sys
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import urlsplit
import zipfile

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
VERSION = '1.0.0-public-rc1'
MAX_UPLOAD = 25 * 1024 * 1024
MAX_REQUEST = 36 * 1024 * 1024
STATIC = {
    '/': 'static/index.html', '/static/index.html': 'static/index.html',
    '/static/public.js': 'static/public.js', '/static/public.css': 'static/public.css',
    '/static/structure.html': 'static/structure.html',
    '/static/structure.css': 'static/structure.css',
    '/static/structure-core.js': 'static/structure-core.js',
    '/static/structure-viewer.js': 'static/structure-viewer.js',
    '/static/vendor/3Dmol-min.js': 'static/vendor/3Dmol-min.js',
    '/static/vendor/3Dmol-LICENSE.txt': 'static/vendor/3Dmol-LICENSE.txt',
    '/static/vendor/3Dmol-min.js.LICENSE.txt': 'static/vendor/3Dmol-min.js.LICENSE.txt',
    '/static/vendor/3Dmol-SOURCE.json': 'static/vendor/3Dmol-SOURCE.json',
}

class PublicStore:
    """Only newly uploaded files and this instance's new results are accessible."""
    def __init__(self, directory):
        self.root = Path(directory).expanduser().resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.files = {}
        self.runs = {}

    def upload(self, filename, content):
        if not isinstance(filename, str) or Path(filename).name != filename or '/' in filename or '\\' in filename:
            raise ValueError('Use a plain filename without a directory.')
        suffix = Path(filename).suffix.lower()
        if suffix not in {'.csv', '.tsv', '.xlsx', '.xls'}:
            raise ValueError('Choose an XLSX, XLS, CSV or TSV table.')
        if not isinstance(content, str):
            raise ValueError('Missing file content.')
        try:
            raw = base64.b64decode(content, validate=True)
        except (ValueError, binascii.Error):
            raise ValueError('Invalid upload encoding.') from None
        if not raw or len(raw) > MAX_UPLOAD:
            raise ValueError('The table must contain data and be at most 25 MB.')
        file_id = secrets.token_hex(16)
        folder = self.root / ('input_' + file_id)
        folder.mkdir(mode=0o700)
        path = folder / filename
        with path.open('xb') as f:
            f.write(raw)
        self.files[file_id] = path
        return {'file_id': file_id, 'filename': filename, **self.inspect(file_id)}

    def _input(self, file_id):
        path = self.files.get(file_id)
        if path is None or path.is_symlink() or not path.is_file():
            raise ValueError('This uploaded table is unavailable in the current session.')
        return path

    def inspect(self, file_id, sheet=None, header=None):
        from proteomics import inspect_file
        return inspect_file(self._input(file_id), sheet, header)

    def run(self, kind, file_id=None, config=None, analysis_id=None, gmt_text=None):
        if kind not in {'analysis', 'quality', 'enrichment'} or not isinstance(config, dict):
            raise ValueError('Choose a supported analysis and a JSON settings object.')
        path = self._input(file_id) if kind != 'enrichment' else None
        if kind == 'enrichment':
            prior = self.runs.get(analysis_id)
            if not prior or prior['kind'] != 'analysis':
                raise ValueError('Run a table analysis before gene-set enrichment.')
            if not isinstance(gmt_text, str) or not gmt_text.strip() or len(gmt_text.encode()) > MAX_UPLOAD:
                raise ValueError('Choose a local GMT file of at most 25 MB.')
        run_id = secrets.token_hex(16)
        out = self.root / ('run_' + run_id)
        out.mkdir(mode=0o700)
        if kind == 'analysis':
            from proteomics import analyze_file
            result = analyze_file(path, config, out)
        elif kind == 'quality':
            from quality import analyze_quality
            result = analyze_quality(path, config, out)
        else:
            from quality import enrich_from_rows
            result = enrich_from_rows(prior['result']['rows'], gmt_text, out, config)
        # No old result directory is scanned or imported.
        self.runs[run_id] = {'kind': kind, 'directory': out, 'result': result}
        return {'run_id': run_id, 'kind': kind, 'summary': result.get('summary', {}),
                'warnings': result.get('warnings', []), 'files': result.get('files', {}),
                'table_preview': result.get('rows', [])[:100]}

    def download(self, run_id, key):
        run = self.runs.get(run_id)
        if not run:
            raise ValueError('This result is unavailable in the current session.')
        if key == 'bundle':
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as archive:
                for name in run['result'].get('files', {}).values():
                    path = self._result_path(run, name)
                    archive.writestr(path.name, path.read_bytes())
            return buf.getvalue(), 'proteomics_' + run_id + '.zip', 'application/zip'
        name = run['result'].get('files', {}).get(key)
        if name is None:
            raise ValueError('Unknown result file.')
        path = self._result_path(run, name)
        return path.read_bytes(), path.name, mimetypes.guess_type(path.name)[0] or 'application/octet-stream'

    @staticmethod
    def _result_path(run, name):
        if not isinstance(name, str) or Path(name).name != name or '/' in name or '\\' in name:
            raise ValueError('Invalid result filename.')
        path = run['directory'] / name
        if path.is_symlink() or not path.is_file() or path.resolve().parent != run['directory'].resolve():
            raise ValueError('Result file is unavailable.')
        return path


def dispatch(store, route, data):
    if not isinstance(data, dict):
        raise ValueError('Request must be a JSON object.')
    if route == '/api/upload':
        return store.upload(data.get('filename'), data.get('content_base64'))
    if route == '/api/inspect':
        return store.inspect(data.get('file_id'), data.get('sheet'), data.get('header'))
    if route in {'/api/analyze', '/api/quality'}:
        return store.run('analysis' if route == '/api/analyze' else 'quality', data.get('file_id'), data.get('config'))
    if route == '/api/enrichment':
        return store.run('enrichment', config=data.get('config'), analysis_id=data.get('analysis_id'), gmt_text=data.get('gmt_text'))
    raise ValueError('This function is not available in this edition.')


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *_):
        pass

    def _local(self):
        origin = 'http://127.0.0.1:' + str(self.server.server_port)
        return (self.headers.get('Host') == origin.removeprefix('http://') and
                self.headers.get('Origin', origin) == origin and
                self.headers.get('Sec-Fetch-Site', 'same-origin') != 'cross-site')

    def _send(self, code, content, mime='application/json', filename=None):
        self.send_response(code)
        self.send_header('Content-Type', mime)
        self.send_header('Content-Length', str(len(content)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('Referrer-Policy', 'no-referrer')
        if filename:
            safe = re.sub(r'[^A-Za-z0-9._-]', '_', filename)
            self.send_header('Content-Disposition', f'attachment; filename="{safe}"')
        self.end_headers()
        self.wfile.write(content)

    def _json(self, code, data):
        self._send(code, json.dumps(data, allow_nan=False).encode())

    def do_GET(self):
        if not self._local():
            return self._json(403, {'error': 'Use this local instance at its printed address.'})
        route = urlsplit(self.path).path
        if route == '/api/session':
            return self._json(200, {'token': self.server.session_token, 'version': VERSION, 'edition': 'public'})
        name = STATIC.get(route)
        if name is None:
            return self._json(404, {'error': 'Not available in this edition.'})
        path = ROOT / name
        if path.is_symlink() or not path.is_file():
            return self._json(404, {'error': 'Public asset unavailable.'})
        return self._send(200, path.read_bytes(), mimetypes.guess_type(path.name)[0] or 'text/plain')

    def do_POST(self):
        if not self._local() or not hmac.compare_digest(self.headers.get('X-Local-Token', ''), self.server.session_token):
            return self._json(403, {'error': 'Reload the local page.'})
        try:
            if self.headers.get('Content-Type', '').split(';')[0] != 'application/json':
                raise ValueError('JSON requests only.')
            size = int(self.headers.get('Content-Length', '-1'))
            if not 0 < size <= MAX_REQUEST:
                raise ValueError('Request is missing or too large.')
            data = json.loads(self.rfile.read(size))
            route = urlsplit(self.path).path
            if route == '/api/download':
                if not isinstance(data, dict):
                    raise ValueError('Request must be a JSON object.')
                content, name, mime = self.server.store.download(data.get('run_id'), data.get('key'))
                return self._send(200, content, mime, name)
            return self._json(200, dispatch(self.server.store, route, data))
        except (ValueError, OSError, TypeError, KeyError) as exc:
            return self._json(400, {'error': str(exc)})
        except Exception:
            return self._json(500, {'error': 'Analysis failed. No successful result was registered.'})


def main(argv=None):
    parser = argparse.ArgumentParser(description='Local proteomics analysis and PDB structure viewer.')
    parser.add_argument('--port', type=int, default=8765)
    parser.add_argument('--data-dir', default='public_data', help='Separate directory for uploaded tables and generated results.')
    parser.add_argument('--check', action='store_true', help='Check installed generic dependencies without starting a service.')
    args = parser.parse_args(argv)
    missing = [x for x in ['numpy', 'pandas', 'scipy', 'openpyxl', 'xlrd', 'plotly', 'lxml'] if importlib.util.find_spec(x) is None]
    if missing:
        print('Missing dependencies: ' + ', '.join(missing), file=sys.stderr)
        return 1
    if args.check:
        print('Generic dependencies found. This is not fresh-install or UI acceptance.')
        return 0
    if not 1 <= args.port <= 65535:
        parser.error('Port must be between 1 and 65535.')
    try:
        server = HTTPServer(('127.0.0.1', args.port), Handler)
    except OSError:
        print('Port unavailable. No existing service was probed, reused, stopped, or opened.', file=sys.stderr)
        return 1
    try:
        server.store = PublicStore(args.data_dir)
        server.session_token = secrets.token_urlsafe(32)
        print(f'Devil\'s Proteomics {VERSION}\nhttp://127.0.0.1:{args.port}/\nData: {server.store.root}', flush=True)
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
