"""Portable offline HTML exports used by proteomics.py."""
import html
import json
from pathlib import Path

from figure_style import math_assets, volcano_formulae

CSS = """
*{box-sizing:border-box}body{font:14px/1.6 Arial,sans-serif;color:#172B4D;background:#F5F7FB;margin:0;padding:28px}body,button,input,select,textarea,pre,code,svg text{font-family:Arial,sans-serif}a{color:#00539B}h1{font-size:28px;color:#012169;margin:0 0 6px}h2{font-size:19px;margin:22px 0 12px}h3{font-size:15px;margin:14px 0 8px}.muted{color:#64748B}.card{background:#fff;border:1px solid #D8E1EE;border-radius:6px;padding:20px;margin:18px 0}.stats{display:grid;grid-template-columns:repeat(4,minmax(110px,1fr));gap:12px}.stat{padding:16px;background:#fff;border:1px solid #D8E1EE;border-radius:6px}.stat strong{font-size:28px;display:block;color:#012169}.notice{padding:10px 14px;background:#EDF3FB;border-left:3px solid #00539B;margin:10px 0}button,input{font:inherit}input{border:1px solid #CBD5E1;border-radius:4px;padding:9px 12px;width:min(500px,100%)}button{cursor:pointer;padding:8px 14px;border-radius:4px;border:1px solid #00539B;background:white;color:#012169}table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:left;vertical-align:top;border-bottom:1px solid #E2E8F0;padding:9px 10px}th{color:#334155;font-weight:600;white-space:nowrap;background:#F8FAFC}td.number{text-align:right;font-variant-numeric:tabular-nums}.scroll{overflow:auto;max-height:600px}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px;color:#334155}.tag{display:inline-block;border-radius:3px;padding:2px 8px;background:#EDF3FB;color:#00539B;font-size:12px}.up{color:#B45309}.down{color:#00539B}.layout{display:grid;grid-template-columns:minmax(0,1fr) 360px;gap:20px}#plot{height:640px}#detail{max-height:750px;overflow:auto}.detailtable td{padding:5px 7px}.chord{width:100%;height:210px}.tooltip{position:fixed;z-index:30;width:350px;max-width:calc(100vw - 24px);max-height:580px;overflow:hidden;pointer-events:none;background:white;border:1px solid #CBD5E1;border-radius:6px;padding:16px;box-shadow:0 8px 20px #01216920;font-size:12px}.tooltip .chord{height:150px}.empty{padding:20px;color:#64748B;text-align:center}.legend{display:flex;gap:18px;color:#64748B;font-size:12px}details{padding:8px 0}summary{cursor:pointer;color:#00539B}.footer{margin-top:20px;color:#64748B;font-size:12px}@media(max-width:900px){body{padding:18px}.layout{grid-template-columns:1fr}.stats{grid-template-columns:1fr 1fr}#plot{height:480px}}@media print{body{background:white;padding:0}.scroll{max-height:none;overflow:visible}.no-print{display:none}}
"""

CSS += """
body{background:#fff;color:#111;max-width:1500px;margin:auto}h1,h2,h3,.muted,.footer,.legend,th,pre,.up,.down{color:#111}
.card,.stat{border-color:#b7b7b7;border-radius:0}.stat strong{color:#111}.notice{background:#fff;border-color:#111}
th{background:#fff;border-top:1px solid #111;border-bottom:1px solid #111}.tag{color:#111;background:#fff;border:1px solid #b7b7b7}
.legend{font-size:14px}.figure-caption{font-size:14px;margin-top:12px}.figure-actions{display:flex;gap:16px;flex-wrap:wrap;margin:12px 0}
.figure-method{border-top:1px solid #b7b7b7;margin-top:18px;padding-top:12px}.figure-method .math-formula{font-size:16px}
button,input{border-radius:0;color:#111}#plot{min-width:0}.math-formula{max-width:100%}.tooltip{border-radius:0;color:#111}
"""

CSS += """
.report-heading{display:flex;justify-content:space-between;align-items:start;gap:24px;flex-wrap:wrap}.report-heading h1{margin:0}.report-source{margin:8px 0 22px;overflow-wrap:anywhere}.report-source span{display:block}.report-kind{border:1px solid #555;padding:6px 12px;font-size:14px}.report-filter{display:flex;gap:16px;align-items:center;flex-wrap:wrap;margin:16px 0}.report-filter input{width:min(560px,100%)}.report-table{max-height:600px;overflow:auto;border-top:1px solid #999}.report-table table{font-size:14px;min-width:800px}.report-table th{position:sticky;top:0;background:white;z-index:1;padding-top:12px;padding-bottom:12px}.report-table td{padding:12px 10px;vertical-align:middle}.report-table tr:hover{background:#f2f5f8}.report-table .select-protein{border:0;background:none;padding:0;text-align:left;color:#00539B;font-weight:700}.report-table .id-label{display:block;font-size:13px;margin-top:3px;color:#111}.report-table .up{color:#A44100}.report-table .down{color:#0072B2}.detail-grid{display:grid;grid-template-columns:minmax(200px,1fr) minmax(200px,1fr);gap:28px}.detail-header{display:flex;justify-content:space-between;align-items:baseline;gap:16px;flex-wrap:wrap}.detail-header h2{margin:0}.detail-header a{font-size:14px}.detail-metrics{display:flex;gap:24px;flex-wrap:wrap;margin:18px 0}.detail-metrics strong{display:block;font-size:20px}.detail-metrics span{font-size:14px}.annotation-section{margin:20px 0}.annotation-section>h3{font-size:18px}.card#detail{max-height:none;overflow:visible;scroll-margin-top:20px}.annotation-preview{height:200px;overflow:hidden}.annotation-preview svg{width:100%;height:100%;font-family:Arial,sans-serif}.tooltip .annotation-preview{height:180px}.input-values{overflow:auto;max-height:300px}.input-values table{table-layout:fixed}.input-values td{overflow-wrap:anywhere}#detail:empty{display:none}@media(max-width:700px){.detail-grid{grid-template-columns:1fr}.report-filter{gap:10px}}@media print{.report-table{max-height:none;overflow:visible}.report-table table{min-width:0}.select-protein{color:#111!important}.pn .pn-canvas{break-inside:avoid}}
"""

def js_json(value):
    return json.dumps(value, ensure_ascii=False, allow_nan=False).replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")

def summary_html(result):
    s = result["summary"]
    items = [("Input rows", "total_rows"), ("Plotted rows", "valid_rows"), ("Significant rows", "significant_rows"), ("Invalid rows", "invalid_rows")]
    return '<div class="stats">' + ''.join(f'<div class="stat"><span class="muted">{label}</span><strong>{s[key]:,}</strong></div>' for label, key in items) + '</div>'

def thresholds(c):
    effect = f'|log{c["x_base"]:g}(FC)| ≥ {c["fc_cutoff"]:g}' if c["fc_cutoff_mode"] == "log" else f'FC ≥ {c["fc_cutoff"]:g} or FC ≤ {1/c["fc_cutoff"]:g}'
    significance = f'{c["significance_metric"]} ≤ {c["significance_cutoff"]:g}' if c["significance_mode"] == "probability" else f'−log{c["y_base"]:g}({c["significance_metric"]}) ≥ {c["significance_cutoff"]:g}'
    return effect + " and " + significance + ". " + {"both": "Both directions", "up": "Up only", "down": "Down only"}[c["direction"]] + ". Inclusive cutoffs."

def document(title, body, scripts):
    return f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)}: Devil\'s Proteomics Public Edition</title><style>{CSS}</style></head><body>{body}<footer style="margin-top:24px;padding-top:16px;border-top:1px solid #999;color:#111;font-size:14px">Made by Zhenyao Sun<br>Department of Chemistry<br>Duke University</footer>{scripts}{math_assets()}</body></html>'

DETAIL_JS = r"""
const esc = v => String(v ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt = v => v == null ? 'N/A' : Number(v).toLocaleString('en-US',{maximumSignificantDigits:5});
function detailHTML(r,compact=false) {
 const name=r.gene||r.protein||'Unnamed protein',abund=Object.entries(r.abundance||{});
 if(compact)return `<h3>${esc(name)}</h3><p>${esc(r.description||'')}</p><table class="detailtable"><tr><td>log FC</td><td>${fmt(r.log_fc)}</td></tr><tr><td>p / q</td><td>${fmt(r.p_value)} / ${fmt(r.q_value)}</td></tr></table><p>${esc(r.function||'No function annotations supplied.')}</p><p>Click to view the supplied record.</p>`;
 return `<div class="detail-header"><h2>${esc(name)}</h2></div><p>${esc(r.protein)}<br>${esc(r.description||'')}</p><div class="detail-metrics"><div><span>log FC</span><strong>${fmt(r.log_fc)}</strong></div><div><span>p value</span><strong>${fmt(r.p_value)}</strong></div><div><span>q value</span><strong>${fmt(r.q_value)}</strong></div><div><span>Status</span><strong>${r.significant?(r.direction==='up'?'Increased':r.direction==='down'?'Decreased':'Significant'):'Not significant'}</strong></div></div><section class="annotation-section"><h3>Supplied function annotations</h3><p>${esc(r.function||'No function annotations supplied.')}</p></section><details><summary>Input abundance (${abund.length} samples)</summary><div class="input-values"><table>${abund.map(([k,v])=>`<tr><td>${esc(k)}</td><td>${fmt(v)}</td></tr>`).join('')}</table></div></details><details><summary>Input values</summary><div class="input-values"><table>${Object.entries(r.source||{}).map(([k,v])=>`<tr><td>${esc(k)}</td><td>${esc(v)}</td></tr>`).join('')}</table></div></details>`;
}
function showProtein(r,scroll=false){
 const detail=document.getElementById('detail');detail.innerHTML=detailHTML(r);
 if(scroll)detail.scrollIntoView({behavior:'smooth',block:'start'});
}
"""

VOLCANO_JS = r"""
const C=DATA.config, valid=DATA.rows.filter(r=>r.valid), plot=document.getElementById('plot'), tip=document.getElementById('tooltip');
const groups=[['Not significant',r=>!r.significant,'#B3B3B3'],['Down',r=>r.significant&&r.direction==='down','#0072B2'],['Up',r=>r.significant&&r.direction==='up','#D55E00'],['No fold change',r=>r.significant&&r.direction==='unchanged','#555555']];
const traces=groups.map(([name,pred,color])=>{const rs=valid.filter(pred);return {type:'scatter',mode:'markers',name:`${name} (n=${rs.length})`,x:rs.map(r=>r.x),y:rs.map(r=>r.y),customdata:rs,marker:{color,size:5,opacity:.85,line:{width:0}},hoverinfo:'none',showlegend:rs.length>0};});
const shapes=[-C.x_threshold,C.x_threshold].filter(Number.isFinite).map(x=>({type:'line',x0:x,x1:x,y0:0,y1:1,yref:'paper',layer:'below',line:{color:'#777777',dash:'dash',width:1}}));if(Number.isFinite(C.y_threshold))shapes.push({type:'line',x0:0,x1:1,xref:'paper',y0:C.y_threshold,y1:C.y_threshold,layer:'below',line:{color:'#777777',dash:'dash',width:1}});
const axis={showline:true,linecolor:'#111',linewidth:1,ticks:'outside',ticklen:5,tickcolor:'#111',showgrid:false,zeroline:false,automargin:true};
const xmax=valid.reduce((m,r)=>Math.max(m,Math.abs(r.x)),Math.max(1,C.x_threshold||0))*1.08;
Plotly.newPlot(plot,traces,{template:'none',margin:{l:80,r:25,b:75,t:60},paper_bgcolor:'#fff',plot_bgcolor:'#fff',font:{family:'Arial',size:13,color:'#111'},legend:{orientation:'h',x:0,y:1.12,font:{size:12,color:'#111'}},xaxis:{...axis,title:{text:`log<sub>${C.x_base}</sub>(treatment / control)`},range:[-xmax,xmax]},yaxis:{...axis,title:{text:`−log<sub>${C.y_base}</sub>(${C.significance_metric})`},rangemode:'tozero'},shapes,hovermode:'closest'},{responsive:true,scrollZoom:true,displaylogo:false,modeBarButtonsToRemove:['sendChartToCloud','sendDataToCloud'],toImageButtonOptions:{format:'svg',filename:'volcano'}});
let pointer={x:120,y:120};document.addEventListener('mousemove',e=>{pointer={x:e.clientX,y:e.clientY};if(!tip.hidden)placeTip();});
function placeTip(){tip.style.left=Math.max(8,Math.min(pointer.x+18,window.innerWidth-tip.offsetWidth-12))+'px';tip.style.top=Math.max(8,Math.min(pointer.y+14,window.innerHeight-tip.offsetHeight-12))+'px';}
plot.on('plotly_hover',e=>{tip.innerHTML=detailHTML(e.points[0].customdata,true);tip.hidden=false;placeTip();});
plot.on('plotly_unhover',()=>tip.hidden=true);plot.on('plotly_click',e=>{showProtein(e.points[0].customdata,true);tip.hidden=true;});
document.getElementById('reset').onclick=()=>Plotly.relayout(plot,{'xaxis.autorange':true,'yaxis.autorange':true});
document.getElementById('search').onkeydown=e=>{if(e.key!=='Enter')return;const q=e.target.value.trim().toLowerCase(),r=valid.find(r=>(r.protein+' '+r.gene).toLowerCase().includes(q));if(!q)return;if(!r){document.getElementById('detail').innerHTML='<p>No matching protein in the plot.</p>';return;}showProtein(r);Plotly.relayout(plot,{'xaxis.range':[r.x-1.5,r.x+1.5],'yaxis.range':[Math.max(0,r.y-2),r.y+2]});};
"""

def write_volcano(path, result):
    from plotly.offline import get_plotlyjs
    c, s = result["config"], result["summary"]
    body = f'''<h1>Volcano plot</h1><p class="muted">Source: {html.escape(s['source_name'])}<br>Sheet: {html.escape(s['sheet'])}</p>
{summary_html(result)}<div class="card"><strong>Significance</strong><p>{html.escape(thresholds(c))}</p>{volcano_formulae(c)}<p>Only the plotted ordinate uses the probability floor. Significance uses the original values.</p><a href="volcano_matlab.zip" download>Download MATLAB volcano</a></div>
<div><div class="card"><input id="search" placeholder="Protein ID or gene; Enter to locate"><button id="reset">Reset view</button><div id="plot"></div><p class="muted">Hover for details. Click to pin.</p></div><aside class="card" id="detail"><div class="empty">Select a protein.</div></aside></div>
<details class="card"><summary>Data notes</summary>{''.join('<p>'+html.escape(x)+'</p>' for x in result['warnings'])}<p>All source rows remain in all_proteins.csv. Duplicate IDs are not merged.</p></details><div id="tooltip" class="tooltip" hidden></div>'''
    scripts = '<script>' + get_plotlyjs() + '</script><script>const DATA=' + js_json({"rows": result["rows"], "config": c, "summary": s}) + ';' + DETAIL_JS + VOLCANO_JS + '</script>'
    Path(path).write_text(document("Volcano plot", body, scripts), encoding="utf-8")

REPORT_JS = r"""
function render(){const query=document.getElementById('filter').value.toLowerCase(),shown=HITS.filter(r=>(r.protein+' '+r.gene+' '+r.description+' '+r.function).toLowerCase().includes(query));document.getElementById('count').textContent=`${shown.length} / ${HITS.length} rows`;document.getElementById('empty').hidden=shown.length>0;document.getElementById('hits').innerHTML=shown.map(r=>`<tr data-row="${r.row_id}"><td>${r.row_id}</td><td><button class="select-protein">${esc(r.gene||r.protein)}</button><span class="id-label">${esc(r.protein)}</span></td><td>${esc(r.description)}</td><td class="number">${fmt(r.x)}</td><td class="number">${fmt(r.p_value)}</td><td class="number">${fmt(r.q_value)}</td><td class="${r.direction}">${r.direction==='up'?'Increased':r.direction==='down'?'Decreased':'Unchanged'}</td></tr>`).join('');}
document.getElementById('filter').addEventListener('input',render);document.getElementById('hits').onclick=e=>{const tr=e.target.closest('tr');if(!tr)return;const r=HITS.find(r=>String(r.row_id)===tr.dataset.row);showProtein(r,true);};render();
"""

def write_report(path, result, manifest):
    c, s = result["config"], result["summary"]
    hits = sorted((row for row in result["rows"] if row["significant"]), key=lambda r: (r["metric_value"], -abs(r["x"]), r["row_id"]))
    method = "Two-sided " + {"welch": "Welch t", "student": "Student equal-variance t", "paired": "paired t", "mannwhitney": "Mann–Whitney U"}[c.get("test_method", "welch")] + " test on log2 abundance; no imputation" if c["mode"] == "replicates" else "Supplied fold changes and significance values"
    body = f'''<header class="report-heading"><h1>Protein analysis report</h1></header><p class="report-source"><span>{html.escape(s['source_name'])}</span><span>Sheet: {html.escape(s['sheet'])}</span></p>
{summary_html(result)}<div class="card"><h2>Results</h2><p><strong>{s['significant_rows']:,} significant rows</strong>: {s['up']:,} up, {s['down']:,} down, {s['unchanged_significant']:,} unchanged.</p><p>{html.escape(thresholds(c))}</p>{volcano_formulae(c)}<p>{s['duplicate_identifier_rows']:,} duplicate identifier rows retained. Counts refer to rows, not unique proteins.</p><div class="figure-actions"><a href="volcano_matlab.zip" download>Download MATLAB volcano</a><a href="manifest.json" download>Analysis settings</a></div></div>
<details class="card"><summary>Methods and data quality</summary><p>{html.escape(method)}.</p><p>q values: {html.escape(s['q_source'])}. BH test population: {s['bh_test_count']:,} valid, nonexcluded p values, including rows with missing FC.</p><p>FC input scale: {html.escape(c['fc_scale'])}. Axis bases: x = {c['x_base']:g}, y = {c['y_base']:g}. Normalization: {html.escape(c['normalization'])}.</p><p>QC exclusions: {s['excluded_rows']:,}. Other invalid rows: {s['invalid_rows']:,}. Full source fields and exclusion reasons: all_proteins.csv.</p>{''.join('<div class="notice">'+html.escape(x)+'</div>' for x in result['warnings'])}</details>'''
    if c["make_volcano"]:
        body += '<section class="card"><h2>Volcano plot</h2><input id="search" placeholder="Protein ID or gene; Enter to locate"><button id="reset">Reset view</button><div id="plot"></div><p class="figure-caption">Each point is one valid protein row. Dashed lines show the selected thresholds. Blue: decreased; vermillion: increased; grey: not significant. Counts and thresholds use the saved analysis. Select a point for its source fields below.</p></section><div id="tooltip" class="tooltip" hidden></div>'
    body += f'''<section class="card"><h2>Significant proteins</h2><div class="report-filter"><input id="filter" class="no-print" aria-label="Filter significant proteins" placeholder="Protein, gene, description or function"><span id="count"></span></div><div class="report-table"><table><thead><tr><th>Row</th><th>Protein / gene</th><th>Description</th><th>log FC</th><th>p value</th><th>q value</th><th>Change</th></tr></thead><tbody id="hits"></tbody></table></div><div id="empty" class="empty" hidden>No significant rows.</div></section>
<section class="card" id="detail"><p>Select a protein to view its supplied abundance, annotations and source fields.</p></section>'''
    scripts = '<script>const REPORT_CONFIG=' + js_json(c) + ';const HITS=' + js_json(hits) + ';' + DETAIL_JS + REPORT_JS + '</script>'
    if c["make_volcano"]:
        from plotly.offline import get_plotlyjs
        scripts += '<script>' + get_plotlyjs() + '</script><script>const DATA=' + js_json({"rows": result["rows"], "config": c}) + ';' + VOLCANO_JS + '</script>'
    Path(path).write_text(document("Protein report", body, scripts), encoding="utf-8")
