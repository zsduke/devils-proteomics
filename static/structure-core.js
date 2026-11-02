/* Local coordinate identity and observed-sequence utilities. No prediction or network code. */
(function (root, factory) {
  'use strict';
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.DevilStructureCore = api;
}(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  const LIMITS = Object.freeze({bytes: 10 * 1024 * 1024, atoms: 100000, models: 20, residues: 10000, annotations: 100, annotationBytes: 1024 * 1024});
  const AA = Object.freeze({ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V'});
  const ELEMENTS = new Set(['H','C','N','O','P','S','SE']);
  const META = new Set(['HEADER','OBSLTE','TITLE','SPLIT','CAVEAT','COMPND','SOURCE','KEYWDS','EXPDTA','NUMMDL','MDLTYP','AUTHOR','REVDAT','SPRSDE','JRNL','REMARK','DBREF','DBREF1','DBREF2','SEQADV','SEQRES','MODRES','HET','HETNAM','HETSYN','FORMUL','HELIX','SHEET','SSBOND','LINK','CISPEP','SITE','CRYST1','ORIGX1','ORIGX2','ORIGX3','SCALE1','SCALE2','SCALE3','MTRIX1','MTRIX2','MTRIX3','ANISOU','SIGATM','SIGUIJ','CONECT','MASTER']);
  const parsedObjects = new WeakSet();
  const own = (o, k) => Object.prototype.hasOwnProperty.call(o, k);
  function need(ok, message) { if (!ok) throw new Error(message); }
  function freeze(value) {
    if (value && typeof value === 'object' && !Object.isFrozen(value)) {
      for (const item of Object.values(value)) freeze(item);
      Object.freeze(value);
    }
    return value;
  }
  function residueKey(modelId, chain, resSeq, iCode) { return JSON.stringify([modelId, chain, resSeq, iCode]); }
  function integer(field, label, line, positive) {
    need(/^\s*-?\d+\s*$/.test(field), `${label} must use decimal PDB numbering (line ${line}).`);
    const n = Number(field.trim());
    need(Number.isSafeInteger(n) && (!positive || n > 0), `Invalid ${label} (line ${line}).`);
    return n;
  }
  function coordinate(field, line) {
    need(/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)$/.test(field.trim()), `Invalid fixed-column coordinate (line ${line}).`);
    const value = Number(field.trim());
    need(Number.isFinite(value), `Non-finite coordinate (line ${line}).`);
    return value;
  }
  function displayResidue(residue) {
    return `Model ${residue.modelId}, chain ${residue.chain || '(blank)'}, ${residue.resSeq}${residue.iCode} ${residue.resName}`;
  }
  function consecutive(a, b) {
    if (a.modelId !== b.modelId || a.chain !== b.chain) return false;
    if (b.resSeq === a.resSeq + 1 && b.iCode === '') return true;
    if (b.resSeq !== a.resSeq || b.iCode === '') return false;
    return (a.iCode === '' && b.iCode === 'A') || (a.iCode !== '' && b.iCode.charCodeAt(0) === a.iCode.charCodeAt(0) + 1);
  }
  function parsePDB(text) {
    need(typeof text === 'string' && text.length > 0, 'Choose a nonempty PDB coordinate file.');
    need(text.length <= LIMITS.bytes && /^[\x09\x0A\x0D\x20-\x7E]*$/.test(text), 'PDB input must be ASCII and at most 10 MiB.');
    need(!/(?:^|\n)\s*(?:data_|loop_|_atom_site\.)/i.test(text), 'mmCIF is not supported. Choose a fixed-column PDB file.');
    const lines = text.replace(/\r\n?/g, '\n').split('\n');
    const models = []; const modelIds = new Set();
    let active = null, explicit = false, implicit = false, ended = false, atoms = 0, residues = 0, ignoredHetatmCount = 0;
    function createModel(serial, isExplicit) {
      const id = String(serial);
      need(!modelIds.has(id), 'Duplicate MODEL serial; model identities must be unique.');
      need(models.length < LIMITS.models, 'At most 20 models are supported.');
      const model = {id, serial, explicit:isExplicit, residues:[], atoms:[], sourceLines:[], residueMap:new Map(), serials:new Set(), completedChains:new Set(), currentChain:null, currentResidue:null, omittedSinceResidue:false};
      modelIds.add(id); models.push(model); return model;
    }
    for (let offset = 0; offset < lines.length; offset++) {
      const line = lines[offset], lineNo = offset + 1;
      if (!line.trim()) continue;
      const record = line.slice(0, 6).trim();
      need(!ended, `Content after END is not supported (line ${lineNo}).`);
      if (record === 'MODEL') {
        need(!implicit && !active, `Nested or mixed implicit/explicit MODEL records (line ${lineNo}).`);
        need(line.slice(0, 6) === 'MODEL ' && line.slice(6, 10).trim() === '' && line.slice(14).trim() === '', `MODEL must use its fixed serial field (line ${lineNo}).`);
        const serial = integer(line.slice(10, 14), 'MODEL serial', lineNo, true);
        explicit = true; active = createModel(serial, true); continue;
      }
      if (record === 'ENDMDL') {
        need(explicit && active && active.atoms.length > 0, `Unmatched or empty ENDMDL (line ${lineNo}).`);
        active = null; continue;
      }
      if (record === 'END') {
        need(!explicit || active === null, 'Explicit models must close with ENDMDL before END.');
        ended = true; continue;
      }
      if (record === 'HETATM') {
        need(!explicit || active !== null, `HETATM outside an explicit MODEL (line ${lineNo}).`);
        if (active && active.currentResidue) active.omittedSinceResidue = true;
        ignoredHetatmCount++; continue;
      }
      if (record === 'TER') {
        if (active && active.currentChain !== null) {
          active.completedChains.add(active.currentChain); active.currentChain = null; active.currentResidue = null;
          active.sourceLines.push('TER');
        }
        continue;
      }
      if (record !== 'ATOM') {
        need(META.has(record), `Unsupported record ${record || '(blank)'} (line ${lineNo}).`);
        continue;
      }
      need(line.slice(0, 6) === 'ATOM  ' && line.length >= 78 && line.length <= 80 && !line.includes('\t'), `ATOM requires fixed PDB columns through the element field (line ${lineNo}).`);
      need(line.slice(16, 17) === ' ', `Alternate conformers are unsupported; select a single conformer before import (line ${lineNo}).`);
      need(line.slice(11, 12) === ' ' && line.slice(20, 21) === ' ' && line.slice(27, 30).trim() === '', `Unsupported ATOM field layout (line ${lineNo}).`);
      need(line.slice(78, 80).trim() === '', `Charged ATOM records are outside this viewer subset (line ${lineNo}).`);
      if (!active) {
        need(!explicit, `ATOM outside an explicit MODEL (line ${lineNo}).`);
        if (!implicit) { active = createModel(1, false); implicit = true; }
        else active = models[0];
      }
      const serial = integer(line.slice(6, 11), 'atom serial', lineNo, true);
      need(!active.serials.has(serial), `Duplicate atom serial within model ${active.id} (line ${lineNo}).`);
      const name = line.slice(12, 16).trim(), resName = line.slice(17, 20).trim();
      const chain = line.slice(21, 22).trim(), resSeq = integer(line.slice(22, 26), 'residue number', lineNo, false), iCode = line.slice(26, 27).trim();
      const element = line.slice(76, 78).trim().toUpperCase();
      need(/^[A-Z0-9'"*]{1,4}$/.test(name) && /^[A-Z0-9]{1,3}$/.test(resName), `Unsupported atom or residue name (line ${lineNo}).`);
      need(/^[A-Za-z0-9]?$/.test(chain) && /^[A-Z]?$/.test(iCode), `Chain must be blank/alphanumeric and insertion code blank/A-Z (line ${lineNo}).`);
      need(ELEMENTS.has(element), `ATOM requires an explicit element H, C, N, O, P, S or SE (line ${lineNo}).`);
      const x = coordinate(line.slice(30, 38), lineNo), y = coordinate(line.slice(38, 46), lineNo), z = coordinate(line.slice(46, 54), lineNo);
      const key = residueKey(active.id, chain, resSeq, iCode);
      let residue = active.residueMap.get(key);
      if (!residue) {
        need(!active.completedChains.has(chain), `Reused chain after a chain boundary or TER is ambiguous (line ${lineNo}).`);
        if (active.currentChain !== null && active.currentChain !== chain) {
          active.completedChains.add(active.currentChain); active.sourceLines.push('TER'); active.currentResidue = null;
        }
        if (active.currentResidue && active.currentChain === chain) {
          const previous = active.currentResidue;
          need(resSeq > previous.resSeq || (resSeq === previous.resSeq && iCode > previous.iCode), `Residue order must increase within each chain (line ${lineNo}).`);
        }
        residues++; need(residues <= LIMITS.residues, 'At most 10,000 observed residues are supported.');
        residue = {key, modelId:active.id, chain, resSeq, iCode, resName, letter:AA[resName] || 'X', coordinateBreakBefore:active.omittedSinceResidue, atomIndices:[], atomNames:new Set()};
        active.residueMap.set(key, residue); active.residues.push(residue);
      } else {
        need(active.currentResidue === residue && active.currentChain === chain, `Repeated noncontiguous residue identity (line ${lineNo}).`);
        need(!active.omittedSinceResidue, `HETATM interleaved within one ATOM residue is unsupported (line ${lineNo}).`);
        need(residue.resName === resName, `Conflicting residue names for one identity (line ${lineNo}).`);
      }
      need(!residue.atomNames.has(name), `Duplicate atom name within a residue (line ${lineNo}).`);
      atoms++; need(atoms <= LIMITS.atoms, 'At most 100,000 ATOM records are supported.');
      const index = active.atoms.length;
      active.atoms.push({index, serial, name, element, x, y, z, residueKey:key, chain, resSeq, iCode, resName});
      residue.atomIndices.push(index); residue.atomNames.add(name); active.serials.add(serial);
      active.currentChain = chain; active.currentResidue = residue; active.omittedSinceResidue = false;
      // Preserve identity/coordinate fixed columns; discard unneeded fields and all source metadata.
      const displayElement = element.length === 2 ? element[0] + element[1].toLowerCase() : element;
      active.sourceLines.push(line.slice(0, 54).padEnd(76, ' ') + displayElement.padStart(2, ' ') + '  ');
    }
    need(!explicit || active === null, 'Explicit MODEL lacks ENDMDL.');
    need(atoms > 0, 'No supported ATOM coordinates were found. HETATM-only files are unsupported.');
    const resultModels = models.map(m => ({id:m.id, serial:m.serial, explicit:m.explicit,
      residues:m.residues.map(({atomNames, ...r}) => r), atoms:m.atoms, pdb:m.sourceLines.join('\n') + '\nEND\n'}));
    const unknown = resultModels.reduce((n, m) => n + m.residues.filter(r => r.letter === 'X').length, 0);
    const warnings = ['Sequence covers observed ATOM residues only. Missing coordinates are not imputed.'];
    if (ignoredHetatmCount) warnings.push(`${ignoredHetatmCount} HETATM records omitted. Ligands, waters and HETATM residues are not displayed.`);
    if (unknown) warnings.push(`${unknown} nonstandard or unknown residues represented as X.`);
    const parsed = freeze({models:resultModels, warnings, ignoredHetatmCount, observedOnly:true, atomCount:atoms, residueCount:residues});
    parsedObjects.add(parsed); return parsed;
  }
  function checked(parsed) { need(parsedObjects.has(parsed), 'Use the result returned by parsePDB.'); return parsed; }
  function modelView(parsed, modelId) {
    const model = checked(parsed).models.find(m => m.id === modelId);
    need(model, 'Choose a model present in this file.');
    return freeze({id:model.id, modelId:model.id, serial:model.serial, explicit:model.explicit, residues:model.residues,
      atoms:model.atoms, atomResidueKeys:model.atoms.map(a => a.residueKey), pdb:model.pdb});
  }
  function selectResidues(parsed, keys) {
    checked(parsed); need(Array.isArray(keys) && keys.length <= LIMITS.residues, 'Selection must be an array of at most 10,000 residue keys.');
    const wanted = new Set(keys);
    need(wanted.size === keys.length && keys.every(k => typeof k === 'string'), 'Selection contains duplicate or invalid residue keys.');
    const ordered = parsed.models.flatMap(m => m.residues).filter(r => wanted.has(r.key));
    need(ordered.length === keys.length, 'Selection contains a residue not present in the loaded structure.');
    return Object.freeze(ordered);
  }
  function sequenceSegments(parsed, keys) {
    const selected = selectResidues(parsed, keys), segments = [];
    const positions = new Map(parsed.models.flatMap(m => m.residues).map((r, i) => [r.key, i]));
    let previous = null;
    for (const r of selected) {
      if (!previous || !consecutive(previous, r) || r.coordinateBreakBefore || positions.get(r.key) !== positions.get(previous.key) + 1) segments.push({modelId:r.modelId, chain:r.chain,
        start:{resSeq:r.resSeq, iCode:r.iCode}, end:{resSeq:r.resSeq, iCode:r.iCode}, sequence:'', residueKeys:[]});
      const segment = segments[segments.length - 1];
      segment.sequence += r.letter; segment.residueKeys.push(r.key); segment.end = {resSeq:r.resSeq, iCode:r.iCode}; previous = r;
    }
    return freeze({observedOnly:true, segments, selectedResidues:selected.length,
      unknownResidues:selected.filter(r => r.letter === 'X').length, omittedCoordinateGaps:true});
  }
  function exportSelection(parsed, keys) {
    const selected = selectResidues(parsed, keys); need(selected.length > 0, 'Select at least one observed residue before export.');
    const result = sequenceSegments(parsed, keys);
    const json = JSON.stringify({schema:'observed-pdb-selection-v1', observedOnly:true,
      residues:selected.map(({atomIndices, ...r}) => r), segments:result.segments,
      note:'Coordinates present in this file only; gaps and missing residues are not filled.'}, null, 2) + '\n';
    const tsv = ['model_id\tchain\tresidue_number\tinsertion_code\tresidue_name\tone_letter\tresidue_key',
      ...selected.map(r => [r.modelId, r.chain, r.resSeq, r.iCode, r.resName, r.letter, r.key].join('\t'))].join('\n') + '\n';
    const fasta = result.segments.map((s, i) => `>segment_${i + 1} model=${s.modelId} chain=${s.chain || '(blank)'} residues=${s.start.resSeq}${s.start.iCode}:${s.end.resSeq}${s.end.iCode} observed_only=true\n${s.sequence.match(/.{1,80}/g).join('\n')}`).join('\n') + '\n';
    return freeze({json, tsv, fasta});
  }
  function plainObject(value) { return value !== null && typeof value === 'object' && !Array.isArray(value) && (Object.getPrototypeOf(value) === Object.prototype || Object.getPrototypeOf(value) === null); }
  function exactFields(value, fields, label) {
    need(plainObject(value) && Object.keys(value).length === fields.length && fields.every(k => own(value, k)), `${label} must contain only ${fields.join(', ')}.`);
  }
  function validateAnnotations(parsed, input) {
    checked(parsed);
    if (typeof input === 'string') {
      need(input.length <= LIMITS.annotationBytes, 'Annotations must be at most 1 MiB.');
      try { input = JSON.parse(input); } catch (_) { throw new Error('Annotations must be valid JSON.'); }
    }
    exactFields(input, ['annotations'], 'Annotation document');
    need(Array.isArray(input.annotations) && input.annotations.length <= LIMITS.annotations, 'Use at most 100 annotations.');
    let total = 0;
    const output = input.annotations.map(item => {
      exactFields(item, ['label','type','residueKeys'], 'Each annotation');
      need(typeof item.label === 'string' && item.label.trim() && item.label.length <= 120 && !/[\x00-\x1F\x7F]/.test(item.label), 'Annotation labels must contain 1–120 printable characters.');
      need(['domain','site','region'].includes(item.type), 'Annotation type must be domain, site or region.');
      const selected = selectResidues(parsed, item.residueKeys); need(selected.length > 0, 'Annotations must identify at least one observed residue.');
      total += selected.length; need(total <= 100000, 'Annotations contain too many residue references.');
      return {label:item.label.trim(), type:item.type, residueKeys:selected.map(r => r.key)};
    });
    return freeze(output);
  }
  return Object.freeze({LIMITS, parsePDB, modelView, selectResidues, sequenceSegments, exportSelection, validateAnnotations, displayResidue});
}));
