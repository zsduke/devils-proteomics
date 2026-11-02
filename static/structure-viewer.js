/* Generic local coordinate display and explicit user selection. No network or prediction. */
(function (root, factory) {
  'use strict';
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.DevilStructureViewer = api;
  if (root.document && root.DevilStructureCore) api.boot(root);
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  const PAGE_SIZE = 100;
  function need(ok, text) { if (!ok) throw new Error(text); }
  function chainName(chain) { return chain === '' ? '(blank)' : chain; }
  function residueNumber(residue) { return String(residue.resSeq) + residue.iCode; }
  function residueLabel(residue) { return 'Model ' + residue.modelId + ' / Chain ' + chainName(residue.chain) + ' / ' + residue.resName + ' ' + residueNumber(residue); }
  function verifyAtomMapping(view, atoms) {
    need(Array.isArray(atoms) && atoms.length === view.atoms.length, 'Rendered atom count differs from the selected model.');
    const keys = [];
    atoms.forEach(function (actual, index) {
      const expected = view.atoms[index];
      const fieldsMatch = actual.index === index && expected.index === index && actual.serial === expected.serial &&
        String(actual.atom || '').trim() === expected.name && String(actual.elem || '').toUpperCase() === expected.element &&
        String(actual.chain || '').trim() === expected.chain && actual.resi === expected.resSeq &&
        String(actual.icode || '').trim() === expected.iCode && actual.resn === expected.resName;
      const coordinatesMatch = ['x', 'y', 'z'].every(axis => Number.isFinite(actual[axis]) && Math.abs(actual[axis] - expected[axis]) <= 0.000001);
      need(fieldsMatch && coordinatesMatch, 'Rendered atom identity differs at index ' + index + '. Selection is disabled.');
      keys.push(expected.residueKey);
    });
    return Object.freeze(keys);
  }
  class SelectionSession {
    constructor(core, parsed) {
      need(parsed.models.length > 0, 'The file contains no supported model.');
      this.core = core; this.parsed = parsed; this.modelId = parsed.models[0].id; this.keys = new Set(); this.annotations = [];
    }
    get view() { return this.core.modelView(this.parsed, this.modelId); }
    chooseModel(id) { this.core.modelView(this.parsed, id); this.modelId = id; this.keys.clear(); return this.view; }
    setSelected(keys) {
      const residues = this.core.selectResidues(this.parsed, keys);
      need(residues.every(r => r.modelId === this.modelId), 'Select residues from the displayed model only.');
      this.keys = new Set(residues.map(r => r.key));
    }
    toggle(key) { const keys = new Set(this.keys); keys.has(key) ? keys.delete(key) : keys.add(key); this.setSelected(Array.from(keys)); }
    clear() { this.keys.clear(); }
    segments() { return this.core.sequenceSegments(this.parsed, Array.from(this.keys)); }
    exports() { return this.core.exportSelection(this.parsed, Array.from(this.keys)); }
    loadAnnotations(input) { this.annotations = this.core.validateAnnotations(this.parsed, input); }
    annotationKeys(annotation) { return this.core.selectResidues(this.parsed, annotation.residueKeys).filter(r => r.modelId === this.modelId).map(r => r.key); }
  }
  function annotationTracks(session) {
    const residues = session.view.residues, chains = Array.from(new Set(residues.map(r => r.chain)));
    return chains.map(function (chain) {
      const observed = residues.filter(r => r.chain === chain), positions = new Map(observed.map((r, i) => [r.key, i]));
      const tracks = session.annotations.map(function (annotation) {
        const keys = session.annotationKeys(annotation).filter(key => positions.has(key));
        const segments = session.core.sequenceSegments(session.parsed, keys).segments;
        return { label: annotation.label, type: annotation.type, runs: segments.map(s => ({ start: positions.get(s.residueKeys[0]), count: s.residueKeys.length, first: residueNumber(s.start), last: residueNumber(s.end) })) };
      }).filter(track => track.runs.length);
      return { chain, observed: observed.length, tracks };
    }).filter(group => group.tracks.length);
  }
  function boot(win) {
    const doc = win.document, core = win.DevilStructureCore, byId = id => doc.getElementById(id);
    let session = null, viewer = null, glModel = null, atomKeys = [], rendered = false;
    let page = 0, fileGeneration = 0, annotationGeneration = 0;
    const exportIds = ['copy-sequence', 'export-fasta', 'export-json', 'export-tsv'];
    function text(id, value) { byId(id).textContent = value; }
    function error(message) { text('error', message); byId('error').hidden = !message; }
    function announce(message) { text('status', message); }
    function guard(action) { return function () { try { error(''); return action.apply(null, arguments); } catch (exc) { error(exc.message || String(exc)); } }; }
    function node(tag, content, className) { const n = doc.createElement(tag); if (content !== undefined) n.textContent = content; if (className) n.className = className; return n; }
    function download(contents, filename, type) {
      const blob = new win.Blob([contents], { type: type || 'text/plain;charset=utf-8' });
      const url = win.URL.createObjectURL(blob), a = node('a'); a.href = url; a.download = filename;
      doc.body.appendChild(a); a.click(); a.remove(); win.setTimeout(() => win.URL.revokeObjectURL(url), 1000);
    }
    function disableInputs(disabled) {
      ['model', 'display', 'chain', 'residue-filter', 'annotation-file', 'annotation-template', 'select-page', 'clear-selection'].forEach(id => { byId(id).disabled = disabled; });
      byId('reset-view').disabled = disabled || !rendered;
    }
    function resetState() {
      annotationGeneration += 1; session = null; rendered = false; glModel = null; atomKeys = []; page = 0;
      if (viewer) { viewer.removeAllModels(); viewer.render(); }
      ['model', 'chain', 'residue-list', 'annotations', 'annotation-map', 'sequence-segments'].forEach(id => byId(id).replaceChildren());
      ['file-summary', 'warnings', 'page-status'].forEach(id => text(id, ''));
      byId('warnings').hidden = true; byId('copy-fallback-label').hidden = true;
      byId('copy-fallback').value = ''; byId('residue-filter').value = ''; byId('annotation-file').value = '';
      text('selection-count', '0 residues selected'); text('viewer-status', 'No coordinates loaded.');
      exportIds.concat(['previous-page', 'next-page']).forEach(id => { byId(id).disabled = true; }); disableInputs(true);
    }
    function filteredResidues() {
      if (!session) return [];
      const query = byId('residue-filter').value.trim().toUpperCase(), chain = byId('chain').value;
      return session.view.residues.filter(r => (chain === '*' || JSON.stringify(r.chain) === chain) && (!query || residueLabel(r).toUpperCase().includes(query) || r.letter === query));
    }
    function pageResidues() { return filteredResidues().slice(page * PAGE_SIZE, (page + 1) * PAGE_SIZE); }
    function drawResidues() {
      const all = filteredResidues(); page = Math.min(page, Math.max(0, Math.ceil(all.length / PAGE_SIZE) - 1));
      const list = byId('residue-list'); list.replaceChildren();
      all.slice(page * PAGE_SIZE, (page + 1) * PAGE_SIZE).forEach(function (residue) {
        const label = node('label', undefined, 'residue-option' + (session.keys.has(residue.key) ? ' selected' : ''));
        const checkbox = node('input'); checkbox.type = 'checkbox'; checkbox.checked = session.keys.has(residue.key);
        checkbox.setAttribute('aria-label', residueLabel(residue));
        checkbox.addEventListener('change', guard(function () { session.toggle(residue.key); updateSelection(); }));
        label.append(checkbox, node('span', 'Chain ' + chainName(residue.chain) + ' / ' + residue.resName + ' ' + residueNumber(residue) + ' / ' + residue.letter)); list.appendChild(label);
      });
      text('page-status', all.length ? 'Showing ' + (page * PAGE_SIZE + 1) + '–' + Math.min((page + 1) * PAGE_SIZE, all.length) + ' of ' + all.length + ' residues' : 'No matching residues.');
      byId('previous-page').disabled = page === 0; byId('next-page').disabled = (page + 1) * PAGE_SIZE >= all.length;
      byId('select-page').disabled = all.length === 0;
    }
    function paintSelection() {
      if (!viewer || !glModel || !rendered) return;
      const base = byId('display').value === 'sticks' ? { stick: { radius: 0.13, color: '#778491' } } :
        { cartoon: { color: '#aab7c5' }, stick: { radius: 0.055, color: '#778491' } };
      glModel.setStyle({}, base);
      const indices = session.view.atoms.filter(a => session.keys.has(a.residueKey)).map(a => a.index);
      if (indices.length) glModel.setStyle({ index: indices }, { cartoon: { color: '#0057b8' }, stick: { color: '#0057b8', radius: 0.2 }, sphere: { color: '#0057b8', scale: 0.22 } });
      viewer.render();
    }
    function updateSelection() {
      drawResidues(); paintSelection();
      const result = session.segments(), target = byId('sequence-segments'); target.replaceChildren();
      text('selection-count', session.keys.size + ' residues selected' + (result.unknownResidues ? ' / ' + result.unknownResidues + ' unknown' : ''));
      result.segments.forEach(function (segment) {
        const box = node('div', undefined, 'segment');
        box.append(node('h4', 'Model ' + segment.modelId + ' / Chain ' + chainName(segment.chain) + ' / ' + residueNumber(segment.start) + '–' + residueNumber(segment.end)), node('p', segment.sequence, 'sequence')); target.appendChild(box);
      });
      if (!result.segments.length) target.appendChild(node('p', 'Select residues in the coordinate view or residue list.'));
      exportIds.forEach(id => { byId(id).disabled = session.keys.size === 0; });
      byId('clear-selection').disabled = session.keys.size === 0; byId('copy-fallback-label').hidden = true;
    }
    function drawAnnotations() {
      const target = byId('annotations'); target.replaceChildren();
      const map = byId('annotation-map'); map.replaceChildren();
      if (!session.annotations.length) { target.appendChild(node('p', 'No annotations loaded.')); return; }
      map.appendChild(node('p', 'Bars follow observed residue order. Numbering gaps and unobserved residues are not filled.'));
      annotationTracks(session).forEach(function (group) {
        const box = node('div', undefined, 'annotation-chain'); box.appendChild(node('h4', 'Chain ' + chainName(group.chain) + ' / ' + group.observed + ' observed residues'));
        group.tracks.forEach(function (track) {
          const row = node('div', undefined, 'annotation-track'); row.appendChild(node('span', track.label));
          const svg = doc.createElementNS('http://www.w3.org/2000/svg', 'svg'); svg.setAttribute('viewBox', '0 0 800 32'); svg.setAttribute('role', 'img');
          svg.setAttribute('aria-label', track.label + ', chain ' + chainName(group.chain) + ', residues ' + track.runs.map(r => r.first + ' to ' + r.last).join('; '));
          const baseline = doc.createElementNS('http://www.w3.org/2000/svg', 'rect'); baseline.setAttribute('x', '0'); baseline.setAttribute('y', '14'); baseline.setAttribute('width', '800'); baseline.setAttribute('height', '3'); baseline.setAttribute('fill', '#b4bdc8'); svg.appendChild(baseline);
          track.runs.forEach(function (run) {
            const rect = doc.createElementNS('http://www.w3.org/2000/svg', 'rect');
            rect.setAttribute('x', String(800 * run.start / group.observed)); rect.setAttribute('y', '6'); rect.setAttribute('width', String(Math.max(1, 800 * run.count / group.observed))); rect.setAttribute('height', '20'); rect.setAttribute('fill', track.type === 'site' ? '#111820' : track.type === 'domain' ? '#003087' : '#0057b8');
            const title = doc.createElementNS('http://www.w3.org/2000/svg', 'title'); title.textContent = run.first + '–' + run.last; rect.appendChild(title); svg.appendChild(rect);
          }); row.appendChild(svg); box.appendChild(row);
        }); map.appendChild(box);
      });
      session.annotations.forEach(function (annotation) {
        const keys = session.annotationKeys(annotation), row = node('div', undefined, 'annotation');
        const detail = node('div'); detail.append(node('div', annotation.label, 'annotation-title'), node('div', annotation.type + ' / ' + keys.length + ' residues in model ' + session.modelId, 'annotation-type'));
        const button = node('button', 'Select region'); button.type = 'button'; button.disabled = keys.length === 0;
        button.addEventListener('click', guard(function () { session.setSelected(keys); updateSelection(); announce('Selected ' + annotation.label + ' in model ' + session.modelId + '.'); }));
        row.append(detail, button); target.appendChild(row);
      });
    }
    function showModel() {
      rendered = false; glModel = null; atomKeys = []; byId('reset-view').disabled = true;
      const view = session.view, chain = byId('chain'); chain.replaceChildren();
      const all = node('option', 'All chains'); all.value = '*'; chain.appendChild(all);
      Array.from(new Set(view.residues.map(r => r.chain))).forEach(function (id) { const option = node('option', 'Chain ' + chainName(id)); option.value = JSON.stringify(id); chain.appendChild(option); });
      page = 0; byId('residue-filter').value = '';
      try {
        need(win.$3Dmol && typeof win.$3Dmol.createViewer === 'function', 'The bundled 3D viewer is unavailable.');
        if (!viewer) viewer = win.$3Dmol.createViewer(byId('viewer'), { backgroundColor: '#ffffff', antialias: true });
        viewer.removeAllModels();
        glModel = viewer.addModel(view.pdb, 'pdb', { keepH: true, doAssembly: false, multimodel: false, computeStruct: false });
        const displayedAtoms = glModel.selectedAtoms({});
        atomKeys = verifyAtomMapping(view, displayedAtoms); rendered = true;
        const displayedModel = glModel, displayedSession = session, displayedModelId = session.modelId, verifiedKeys = atomKeys;
        glModel.setClickable({}, true, guard(function (atom) {
          if (!rendered || glModel !== displayedModel || session !== displayedSession || session.modelId !== displayedModelId) return;
          need(Number.isInteger(atom.index) && verifiedKeys[atom.index] && displayedAtoms[atom.index] === atom, 'The clicked atom has no verified residue identity.');
          session.toggle(verifiedKeys[atom.index]); updateSelection();
        }));
        paintSelection(); viewer.zoomTo(); viewer.resize(); viewer.render(); byId('reset-view').disabled = false;
        text('viewer-status', 'Model ' + session.modelId + ' / ' + view.atoms.length + ' verified atoms / ' + view.residues.length + ' observed residues');
      } catch (exc) {
        rendered = false; glModel = null; atomKeys = [];
        if (viewer) { viewer.removeAllModels(); viewer.render(); }
        text('viewer-status', '3D display unavailable. The residue list and sequence exports remain available.'); error(exc.message || String(exc));
      }
      drawAnnotations(); updateSelection();
    }
    async function loadFile() {
      const generation = ++fileGeneration, file = byId('pdb-file').files[0]; resetState(); error('');
      if (!file) { announce('Choose a PDB file.'); return; }
      try {
        need(/\.(pdb|ent)$/i.test(file.name), 'Choose a .pdb or .ent file. mmCIF is not supported.');
        need(file.size > 0 && file.size <= 10 * 1024 * 1024, 'The PDB file must be nonempty and at most 10 MB.');
        const contents = await file.text(); if (generation !== fileGeneration) return;
        const parsed = core.parsePDB(contents); session = new SelectionSession(core, parsed);
        parsed.models.forEach(function (model) { const option = node('option', 'Model ' + model.id + (model.explicit ? '' : ' (implicit)')); option.value = model.id; byId('model').appendChild(option); });
        text('file-summary', file.name + ' / ' + parsed.models.length + ' model' + (parsed.models.length === 1 ? '' : 's'));
        text('warnings', parsed.warnings.join(' ')); byId('warnings').hidden = parsed.warnings.length === 0;
        disableInputs(false); showModel(); announce('Coordinates loaded. Select residues to define a region.');
      } catch (exc) { if (generation !== fileGeneration) return; resetState(); error(exc.message || String(exc)); announce('The file was not loaded.'); }
    }
    byId('pdb-file').addEventListener('change', loadFile);
    byId('model').addEventListener('change', guard(function () { session.chooseModel(byId('model').value); showModel(); announce('Model changed. Previous residue selection cleared.'); }));
    byId('display').addEventListener('change', guard(paintSelection));
    byId('reset-view').addEventListener('click', guard(function () { if (viewer && rendered) { viewer.zoomTo(); viewer.resize(); viewer.render(); } }));
    ['chain', 'residue-filter'].forEach(id => byId(id).addEventListener(id === 'chain' ? 'change' : 'input', guard(function () { page = 0; drawResidues(); })));
    byId('previous-page').addEventListener('click', guard(function () { page -= 1; drawResidues(); }));
    byId('next-page').addEventListener('click', guard(function () { page += 1; drawResidues(); }));
    byId('select-page').addEventListener('click', guard(function () { session.setSelected(Array.from(new Set([...session.keys, ...pageResidues().map(r => r.key)]))); updateSelection(); }));
    byId('clear-selection').addEventListener('click', guard(function () { session.clear(); updateSelection(); }));
    byId('annotation-file').addEventListener('change', async function () {
      const generation = ++annotationGeneration, current = session, file = byId('annotation-file').files[0]; error('');
      if (!file || !current) return;
      try {
        need(file.size > 0 && file.size <= 1024 * 1024, 'Annotation JSON must be nonempty and at most 1 MB.');
        const contents = await file.text(); if (current !== session || generation !== annotationGeneration) return;
        current.loadAnnotations(contents); drawAnnotations(); announce(current.annotations.length + ' annotations loaded.');
      } catch (exc) { if (current === session && generation === annotationGeneration) error(exc.message || String(exc)); }
    });
    byId('annotation-template').addEventListener('click', guard(function () {
      const keys = session.keys.size ? Array.from(session.keys) : session.view.residues.slice(0, 1).map(r => r.key);
      download(JSON.stringify({ annotations: [{ label: 'Selected region', type: 'region', residueKeys: keys }] }, null, 2) + '\n', 'annotations-template.json', 'application/json');
    }));
    [['export-fasta', 'fasta', 'fasta', 'text/plain'], ['export-json', 'json', 'json', 'application/json'], ['export-tsv', 'tsv', 'tsv', 'text/tab-separated-values']].forEach(function (entry) {
      byId(entry[0]).addEventListener('click', guard(function () { download(session.exports()[entry[1]], 'structure-selection.' + entry[2], entry[3] + ';charset=utf-8'); }));
    });
    byId('copy-sequence').addEventListener('click', async function () {
      try {
        const fasta = session.exports().fasta;
        if (win.navigator.clipboard && typeof win.navigator.clipboard.writeText === 'function') {
          try { await win.navigator.clipboard.writeText(fasta); announce('FASTA copied.'); return; } catch (_) { /* Local manual-copy fallback. */ }
        }
        byId('copy-fallback').value = fasta; byId('copy-fallback-label').hidden = false; byId('copy-fallback').focus(); byId('copy-fallback').select();
        announce('FASTA selected below. Use Copy to place it on the clipboard.');
      } catch (exc) { error(exc.message || String(exc)); }
    });
    if (win.location.protocol === 'file:') {
      const link = byId('proteomics-link'); link.removeAttribute('href'); link.setAttribute('aria-disabled', 'true'); link.title = 'Start the local public launcher to use Proteomics.';
    }
    if (typeof win.ResizeObserver === 'function') {
      new win.ResizeObserver(function () { if (viewer && rendered) { viewer.resize(); viewer.render(); } }).observe(byId('viewer'));
    } else win.addEventListener('resize', function () { if (viewer && rendered) { viewer.resize(); viewer.render(); } });
    return Object.freeze({ getSession: () => session, isRendered: () => rendered });
  }
  return Object.freeze({ verifyAtomMapping, SelectionSession, annotationTracks, chainName, residueNumber, residueLabel, boot });
});
