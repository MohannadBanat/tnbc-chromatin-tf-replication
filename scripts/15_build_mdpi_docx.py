"""Build the Biomolecules (MDPI) submission .docx from the manuscript markdown export and the journal template.

Usage: python 15_build_mdpi_docx.py <manuscript.md> <template_unzipped_dir> <figures_dir> <out.docx>

- Section order changed to MDPI's: 1 Introduction, 2 Materials and Methods, 3 Results, 4 Discussion, 5 Conclusions
  (the working draft used Results-before-Methods). Headings and "Section x.y" cross-references are renumbered.
- Figures and tables renumbered in order of first citation; references renumbered in order of first citation
  (MDPI rule), counting citations in text, tables and legends.
- Template styles are used throughout (MDPI12title, MDPI31text, MDPI51figurecaption, MDPI81references, ...).
"""
import sys, re, os, shutil, html, zipfile
from PIL import Image

MD, TPL, FIGS, OUT = sys.argv[1:5]
md = open(MD, encoding='utf-8').read()
md = md.replace('\\[', '[').replace('\\]', ']').replace('\\*', '\u0000STAR\u0000').replace('\\_', '_').replace('\\~', '~')

# ---------------------------------------------------------------- parse markdown into blocks
blocks, para = [], []
lines = md.split('\n'); i = 0
while i < len(lines):
    ln = lines[i]
    if ln.startswith('|'):
        rows = []
        while i < len(lines) and lines[i].startswith('|'):
            rows.append(lines[i]); i += 1
        cells = [[c.strip() for c in r.strip().strip('|').split('|')] for r in rows]
        cells = [r for r in cells if not all(re.fullmatch(r':?-{3,}:?', c) for c in r)]
        blocks.append(('table', cells)); continue
    m = re.match(r'^(#{1,3}) (.*)', ln)
    if m: blocks.append((f'h{len(m.group(1))}', m.group(2))); i += 1; continue
    m = re.match(r'^!\[([^\]]*)\]\(([^)]+)\)', ln) or re.match(r'^&#91;image: ([^\]]*)\]', ln)
    if m: blocks.append(('img', m.group(1))); i += 1; continue
    m = re.match(r'^- (.*)', ln)
    if m: blocks.append(('bullet', m.group(1))); i += 1; continue
    m = re.match(r'^(\d+)[.)] (.*)', ln)
    if m: blocks.append(('ol', (int(m.group(1)), m.group(2)))); i += 1; continue
    if ln.strip(): blocks.append(('p', ln.strip()))
    i += 1

# ---------------------------------------------------------------- split into sections
def h2num(t):
    m = re.match(r'(\d+)\. ', t); return int(m.group(1)) if m else None
sections, cur = {}, 'front'; sections[cur] = []
for b in blocks:
    if b[0] == 'h2':
        n = h2num(b[1]); cur = n if n else b[1]; sections[cur] = [b]
    else: sections[cur].append(b)
concl = sections[5]
backmatter = [b for b in concl if b[0] == 'p' and re.match(r'\*\*[^*]+:\*\*', b[1])]
sections[5] = [b for b in concl if b not in backmatter]
refs = {n: t for k, (n, t) in [(b[0], b[1]) for b in sections['References'] if b[0] == 'ol']}
refnote = [b for b in sections['References'] if b[0] == 'p']
title = next(b[1] for b in sections['front'] if b[0] == 'h1')
abstract = next(b[1] for b in sections['Abstract'] if b[0] == 'p' and not b[1].startswith('**Keywords'))
keywords = next(b[1] for b in sections['Abstract'] if b[0] == 'p' and b[1].startswith('**Keywords')).split(':**', 1)[1].strip()

SEC = {1: 1, 4: 2, 2: 3, 3: 4, 5: 5}
def renum_heading(t):
    return re.sub(r'^(\d+)(\.)', lambda m: f'{SEC[int(m.group(1))]}{m.group(2)}', t)
def renum_sections(t):
    return re.sub(r'\bSection (\d+)\.(\d+)', lambda m: f'Section {SEC[int(m.group(1))]}.{m.group(2)}', t)
body = []
for n in [1, 4, 2, 3, 5]:
    for b in sections[n]:
        if b[0] in ('h2', 'h3'): body.append((b[0], renum_heading(b[1])))
        else: body.append(b)
# data-sources table (Methods) has no caption: cite and caption it
out = []
for j, b in enumerate(body):
    if b[0] == 'table' and (j == 0 or not (body[j - 1][0] == 'p' and body[j - 1][1].startswith('**Table'))):
        out.append(('p', 'The public datasets analysed are listed in Table D.'))
        out.append(('p', '**Table D.** Public datasets analysed and the hypotheses they address.'))
    out.append(b)
body = out

def textof(b):
    if b[0] == 'table': return ' '.join(' '.join(r) for r in b[1])
    if b[0] == 'ol': return b[1][1]
    return b[1] if isinstance(b[1], str) else ''

# ---------------------------------------------------------------- figure / table renumbering by first citation
def first_order(pattern, captions_prefix):
    order = []
    for b in body:
        t = textof(b)
        if b[0] == 'p' and t.startswith(captions_prefix): continue
        for m in re.finditer(pattern, t):
            k = m.group(1)
            if k not in order: order.append(k)
    return order
fig_order = first_order(r'\bFigure (\d+)', '**Figure')
tab_order = first_order(r'\bTable (\d+|D)\b', '**Table')
FIG = {k: str(i + 1) for i, k in enumerate(fig_order)}
TAB = {k: str(i + 1) for i, k in enumerate(tab_order)}
print('figure map', FIG, '| table map', TAB)
def renum_ft(t):
    t = re.sub(r'\bFigure (\d+)', lambda m: f'Figure \u0001{FIG[m.group(1)]}', t)
    t = re.sub(r'\bTable (\d+|D)\b', lambda m: f'Table \u0001{TAB[m.group(1)]}', t)
    return t.replace('\u0001', '')
body = [(k, renum_ft(renum_sections(v))) if isinstance(v, str) else
        (k, [[renum_ft(c) for c in r] for r in v]) if k == 'table' else (k, v) for k, v in body]
backmatter = [(k, renum_ft(renum_sections(v))) for k, v in backmatter]
abstract = renum_sections(abstract)
# (image blocks' alt text "Figure N" was renumbered by renum_ft above)

# ---------------------------------------------------------------- reference renumbering by first citation
CIT = re.compile(r'\[(\d+(?:\s*[,–-]\s*\d+)*)\]')
def expand(s):
    out = []
    for part in re.split(r'\s*,\s*', s):
        if re.search(r'[–-]', part):
            a, b2 = [int(x) for x in re.split(r'[–-]', part)]; out += list(range(a, b2 + 1))
        else: out.append(int(part))
    return out
order = []
for t in [abstract] + [textof(b) for b in body] + [b[1] for b in backmatter]:
    for m in CIT.finditer(t):
        for n in expand(m.group(1)):
            if n not in order: order.append(n)
REF = {old: i + 1 for i, old in enumerate(order)}
uncited = sorted(set(refs) - set(order)); missing = sorted(set(order) - set(refs))
print(f'references cited: {len(order)} of {len(refs)}; uncited (dropped): {uncited}; cited but missing: {missing}')
def compress(nums):
    nums = sorted(set(nums)); out, s = [], 0
    while s < len(nums):
        e = s
        while e + 1 < len(nums) and nums[e + 1] == nums[e] + 1: e += 1
        out.append(f'{nums[s]}–{nums[e]}' if e - s >= 2 else ','.join(str(x) for x in nums[s:e + 1])); s = e + 1
    return ','.join(out)
def renum_cit(t): return CIT.sub(lambda m: '[' + compress([REF[n] for n in expand(m.group(1))]) + ']', t)
body = [(k, renum_cit(v)) if isinstance(v, str) else (k, [[renum_cit(c) for c in r] for r in v]) if k == 'table' else (k, v) for k, v in body]
backmatter = [(k, renum_cit(v)) for k, v in backmatter]
new_refs = [refs[old] for old in order]

# ---------------------------------------------------------------- XML helpers
def esc(s): return html.escape(s, quote=False).replace('\u0000STAR\u0000', '*')
def runs(t, base=''):
    """markdown inline (**bold**, *italic*, `code`) -> w:r elements"""
    out = []
    for tok in re.split(r'(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)', t):
        if not tok: continue
        rpr = base
        if tok.startswith('**'): tok = tok[2:-2]; rpr += '<w:b/>'
        elif tok.startswith('*'): tok = tok[1:-1]; rpr += '<w:i/>'
        elif tok.startswith('`'): tok = tok[1:-1]
        # nested italic inside bold, e.g. **Table 1.** keeps as is
        out.append(f'<w:r>{"<w:rPr>" + rpr + "</w:rPr>" if rpr else ""}<w:t xml:space="preserve">{esc(tok)}</w:t></w:r>')
    return ''.join(out)
def P(style, t, ppr_extra='', base=''):
    return f'<w:p><w:pPr><w:pStyle w:val="{style}"/>{ppr_extra}</w:pPr>{runs(t, base)}</w:p>'

# images
TEXT_W_TW = 7857; EMU = TEXT_W_TW * 635
media = os.path.join(TPL, 'word', 'media'); rels_p = os.path.join(TPL, 'word', '_rels', 'document.xml.rels')
rels = open(rels_p).read()
FIGFILE = {'1': 'Figure1_PDX_retest.png', '2': 'Figure2_dependency.png', '3': 'Figure3_drug_replication.png',
           '4': 'Figure4_activity_reproducibility.png', '5': 'Figure5_crosslab_chromatin.png'}  # draft numbering
INV = {v: k for k, v in FIG.items()}
rid_n, docpr = 900, 5000
def figure_xml(newnum):
    global rid_n, docpr
    src = os.path.join(FIGS, FIGFILE[INV[newnum]]); rid_n += 1; docpr += 1
    name = f'fig{newnum}.png'; shutil.copy(src, os.path.join(media, name))
    global rels
    rels = rels.replace('</Relationships>', f'<Relationship Id="rId{rid_n}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/{name}"/></Relationships>')
    w, h = Image.open(src).size; small = w < 1500   # narrow single-panel figure: 60% of the text column
    FULL = (2608 + TEXT_W_TW) * 635                    # full page width (template allows multi-panel figures across it)
    wide = w / 300 >= 6.0                               # figures drawn >= 6 in wide: centre across the full width at ~true size
    cx = int(min(FULL, w / 300 * 914400)) if wide else int(EMU * (0.6 if small else 1.0)); cy = int(cx * h / w)
    ind = '<w:ind w:left="0"/>' if wide else '<w:ind w:left="2608"/>'
    return (f'<w:p><w:pPr><w:pStyle w:val="MDPI52figure"/>{ind}</w:pPr><w:r><w:rPr><w:noProof/></w:rPr><w:drawing>'
            f'<wp:inline distT="0" distB="0" distL="0" distR="0"><wp:extent cx="{cx}" cy="{cy}"/><wp:effectExtent l="0" t="0" r="0" b="0"/>'
            f'<wp:docPr id="{docpr}" name="Figure {newnum}"/><wp:cNvGraphicFramePr><a:graphicFrameLocks xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" noChangeAspect="1"/></wp:cNvGraphicFramePr>'
            f'<a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
            f'<pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture"><pic:nvPicPr><pic:cNvPr id="{docpr}" name="{name}"/><pic:cNvPicPr/></pic:nvPicPr>'
            f'<pic:blipFill><a:blip r:embed="rId{rid_n}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill><pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{cx}" cy="{cy}"/></a:xfrm>'
            f'<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>')

def table_xml(rows):
    nc = len(rows[0])
    col = lambda c: [re.sub(r'[*`]', '', r[c]) if c < len(r) else '' for r in rows]
    # column width: never narrower than the column's longest word (~95 twips per character at 9 pt), extra space by mean text length
    minw = [max(max((len(w) for w in x.split()), default=1) for x in col(c)) * 95 + 170 for c in range(nc)]
    extra = max(TEXT_W_TW - sum(minw), 0); mean = [sum(len(x) for x in col(c)) / len(rows) + 1 for c in range(nc)]
    ws = [int(minw[c] + extra * mean[c] / sum(mean)) if sum(minw) <= TEXT_W_TW else int(TEXT_W_TW * minw[c] / sum(minw)) for c in range(nc)]
    ws[-1] = TEXT_W_TW - sum(ws[:-1])
    x = (f'<w:tbl><w:tblPr><w:tblW w:w="{TEXT_W_TW}" w:type="dxa"/><w:tblInd w:w="2608" w:type="dxa"/><w:tblBorders>'
         '<w:top w:val="single" w:sz="8" w:space="0" w:color="auto"/><w:bottom w:val="single" w:sz="8" w:space="0" w:color="auto"/></w:tblBorders>'
         '<w:tblLayout w:type="fixed"/><w:tblCellMar><w:left w:w="57" w:type="dxa"/><w:right w:w="57" w:type="dxa"/></w:tblCellMar></w:tblPr><w:tblGrid>'
         + ''.join(f'<w:gridCol w:w="{w}"/>' for w in ws) + '</w:tblGrid>')
    for ri, r in enumerate(rows):
        x += '<w:tr>'
        for c in range(nc):
            border = '<w:tcBorders><w:bottom w:val="single" w:sz="4" w:space="0" w:color="auto"/></w:tcBorders>' if ri == 0 else ''
            txt = r[c] if c < len(r) else ''
            x += (f'<w:tc><w:tcPr><w:tcW w:w="{ws[c]}" w:type="dxa"/>{border}<w:vAlign w:val="center"/></w:tcPr>'
                  + P('MDPI42tablebody', txt, '<w:spacing w:line="240" w:lineRule="auto"/>', '<w:b/>' if ri == 0 else '') + '</w:tc>')
        x += '</w:tr>'
    return x + '</w:tbl>'

# ---------------------------------------------------------------- body XML
xml = []
prev = None
for k, v in body:
    if k == 'h2': xml.append(P('MDPI21heading1', v))
    elif k == 'h3': xml.append(P('MDPI22heading2', v))
    elif k == 'img': xml.append(figure_xml(re.search(r'(\d+)', v).group(1)))
    elif k == 'table': xml.append(table_xml(v))
    elif k == 'bullet': xml.append(P('MDPI38bullet', v))
    elif k == 'p':
        if v.startswith('**Figure'): xml.append(P('MDPI51figurecaption', v))
        elif v.startswith('**Table'): xml.append(P('MDPI41tablecaption', v))
        elif prev and prev[0] == 'table': xml.append(P('MDPI43tablefooter', v))
        else: xml.append(P('MDPI31text', v))
    prev = (k, v)
for k, v in backmatter:
    xml.append(P('MDPI62backmatter', v, '<w:spacing w:before="240"/>' if v.startswith('**Supplementary') else ''))
xml.append(P('MDPI21heading1', 'References'))
for p in refnote: xml.append(P('MDPI31text', p[1]))
for i, t in enumerate(new_refs, 1):
    xml.append(P('MDPI81references', f'{i}.\t{t}', '<w:numPr><w:ilvl w:val="0"/><w:numId w:val="0"/></w:numPr><w:ind w:left="425" w:hanging="425"/>'))
body_xml = ''.join(xml)

# ---------------------------------------------------------------- front matter into the template
doc = open(os.path.join(TPL, 'word', 'document.xml'), encoding='utf-8').read()
def replace_para(doc, style, new, nth=0):
    ms = [m for m in re.finditer(r'<w:p [^>]*>(?:(?!</w:p>).)*?w:pStyle w:val="%s".*?</w:p>' % style, doc, flags=re.S)]
    m = ms[nth]; return doc[:m.start()] + new + doc[m.end():]
SUP = lambda s: f'<w:r><w:rPr><w:vertAlign w:val="superscript"/></w:rPr><w:t xml:space="preserve">{esc(s)}</w:t></w:r>'
T = lambda s, b=False: f'<w:r>{"<w:rPr><w:b/></w:rPr>" if b else ""}<w:t xml:space="preserve">{esc(s)}</w:t></w:r>'
doc = replace_para(doc, 'MDPI11articletype', P('MDPI11articletype', 'Article'))
doc = replace_para(doc, 'MDPI12title', P('MDPI12title', title))
# authors and affiliations as in the authors' companion MDPI manuscript (Cancers), same order
doc = replace_para(doc, 'MDPI13authornames',
    '<w:p><w:pPr><w:pStyle w:val="MDPI13authornames"/></w:pPr>' + T('Areen Arabiat ') + SUP('1') + T(', Alia Abu Kiwan ') + SUP('2')
    + T(', Mohannad Banat ') + SUP('3') + T(', Manal A. Abbas ') + SUP('2,*') + T(' and Hamza Abu Owida ') + SUP('3,*') + '</w:p>')
AFF = '<w:p><w:pPr><w:pStyle w:val="MDPI16affiliation"/></w:pPr>{}<w:r><w:tab/><w:t xml:space="preserve">{}</w:t></w:r></w:p>'
doc = replace_para(doc, 'MDPI16affiliation', AFF.format(SUP('1'), esc('Department of Communications and Computer Engineering, Faculty of Engineering, Al-Ahliyya Amman University, Amman, Jordan')), 0)
doc = replace_para(doc, 'MDPI16affiliation', AFF.format(SUP('2'), esc('Department of Medical Laboratory Sciences, Faculty of Allied Medical Sciences, Al-Ahliyya Amman University, Amman, Jordan'))
                   + AFF.format(SUP('3'), esc('Department of Medical Engineering, Faculty of Engineering, Al-Ahliyya Amman University, Amman, Jordan')), 0)
doc = replace_para(doc, 'MDPI16affiliation', AFF.format('<w:r><w:rPr><w:b/></w:rPr><w:t>*</w:t></w:r>', esc('Correspondence: m.abbas@ammanu.edu.jo (M.A.); h.abuowida@ammanu.edu.jo (H.A.O.)')), 0)
doc = replace_para(doc, 'MDPI17abstract', P('MDPI17abstract', abstract, '<w:spacing w:before="0" w:after="0"/>'), 1)
doc = replace_para(doc, 'MDPI18keywords', '<w:p><w:pPr><w:pStyle w:val="MDPI18keywords"/></w:pPr>' + T('Keywords: ', True) + T(keywords) + '</w:p>')
# (replaced paragraphs are written as plain <w:p>, so each replace_para above targets the next unreplaced one)
# replace template body: from "0. How to Use This Template" heading to just before the publisher's note
start = doc.rfind('<w:p ', 0, doc.find('0. How to Use This Template'))
notes = doc.find('MDPI63notes'); end = doc.rfind('<w:p ', 0, notes)
doc = doc[:start] + body_xml + doc[end:]
# the template's section carries an RTL flag (<w:bidi/>) that Word ignores for English text but other renderers apply to tables
# the copyright/history box is a table floating at a fixed page position; at the template's y = 12541 it overruns the bottom margin
# once filled, and some renderers then push it (and the abstract) to page 2. Raise it slightly.
doc = doc.replace('w:tblpY="12541"', 'w:tblpY="11800"')
doc = doc.replace('<w:bidi/></w:sectPr>', '</w:sectPr>').replace('<w:bidi/><w:docGrid', '<w:docGrid')
open(os.path.join(TPL, 'word', 'document.xml'), 'w', encoding='utf-8').write(doc)
open(rels_p, 'w').write(rels)
ct = os.path.join(TPL, '[Content_Types].xml')
ctx = open(ct).read().replace('wordprocessingml.template.main+xml', 'wordprocessingml.document.main+xml'); open(ct, 'w').write(ctx)
# drop the template-only keyboard-customisation part (not valid in a .docx)
cu = os.path.join(TPL, 'word', 'customizations.xml')
if os.path.exists(cu): os.remove(cu)
ctx = open(ct).read(); ctx = re.sub(r'<Override PartName="/word/customizations.xml"[^>]*/>', '', ctx); open(ct, 'w').write(ctx)
rx = open(rels_p).read(); rx = re.sub(r'<Relationship [^>]*Target="customizations.xml"[^>]*/>', '', rx); open(rels_p, 'w').write(rx)
# zip
if os.path.exists(OUT): os.remove(OUT)
with zipfile.ZipFile(OUT, 'w', zipfile.ZIP_DEFLATED) as z:
    for rootd, _, files in os.walk(TPL):
        for f in files:
            p = os.path.join(rootd, f); z.write(p, os.path.relpath(p, TPL))
print('written', OUT)
