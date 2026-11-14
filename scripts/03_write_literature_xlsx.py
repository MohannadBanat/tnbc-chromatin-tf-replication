"""Write literature_log.xlsx: README, Search log (with COUNTIFS links to Screening), Screening, Key papers."""
import json, sys, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
lit = sys.argv[1]
log = pd.read_csv(f'{lit}/literature_log_draft.csv', dtype=str).fillna('')
counts = json.load(open(f'{lit}/search_counts_2026-11-05.json'))
plan = open('LITERATURE_SEARCH_PLAN.md').read()
import re
strings = dict(re.findall(r'\*\*(S\d+[ab]?)[^*]*\*\*[^\n]*\n(?:PubMed|Also): `([^`]+)`', plan))
strings['S8b'] = re.search(r'Also: `([^`]+)`', plan).group(1)
F = 'Arial'; H = Font(name=F, bold=True, color='FFFFFF'); HF = PatternFill('solid', fgColor='1F3864')
N = Font(name=F, size=10); B = Font(name=F, size=10, bold=True)
def header(ws, cols, widths):
    for i, (c, w) in enumerate(zip(cols, widths), 1):
        x = ws.cell(1, i, c); x.font = H; x.fill = HF; x.alignment = Alignment(wrap_text=True, vertical='top')
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = 'A2'
wb = Workbook()
# README
ws = wb.active; ws.title = 'README'
lines = [
 ('Literature log – TNBC TF-dependency replication study (Biomolecules, MDPI)', B),
 ('Created 2026-11-05 by Claude from PubMed (E-utilities), OpenAlex and Crossref/DataCite.', N),
 ('', N),
 ('What each sheet holds', B),
 ('Search log – the 11 PubMed searches (S1–S10, S8 split into a/b): exact string, date, hits; kept counts are live COUNTIFS on Screening.', N),
 ('Screening – every record (984 PubMed + 26 from citation chasing / methods) with keep/drop, priority, reason and target paper section.', N),
 ('Key papers – the kept records only, sorted by priority (A = must cite, B = likely cite, C = optional).', N),
 ('', N),
 ('How screening was done', B),
 ('Titles screened for all 984 PubMed records; abstracts read only where the title was borderline. Decisions are a first pass and need a human check.', N),
 ('Citation chasing: Bahl et al. 2026 has 0 citations in OpenAlex and Semantic Scholar (checked 2026-11-05); Google Scholar not checked (no API) – check by hand.', N),
 ('All 98 kept DOIs resolved on 2026-11-05 (97 Crossref, 1 DataCite: the Code Ocean capsule). 1 kept record (PMID 41726490, AMIA proceedings) has no DOI.', N),
 ('', N),
 ('How to edit', B),
 ('Change Decision (keep/drop), Priority (A/B/C), Reason and Paper section in the Screening sheet only; the kept counts in Search log update automatically.', N),
 ('Add new records as new rows at the bottom of Screening, using the same columns.', N),
]
for i, (t, f) in enumerate(lines, 1):
    c = ws.cell(i, 1, t); c.font = f; c.alignment = Alignment(wrap_text=True)
ws.column_dimensions['A'].width = 120
# Screening
cols = ['ID','Source','Search','Citation','Title','Year','DOI','PMID','Decision','Priority','Reason / role','Paper section','Screen level']
sc = wb.create_sheet('Screening'); header(sc, cols, [16,14,12,30,70,7,30,11,9,8,50,28,18])
pri = {'A':0,'B':1,'C':2,'':3}
log['_o'] = log.Decision.map({'keep':0,'drop':1}); log['_p'] = log.Priority.map(pri)
log = log.sort_values(['_o','_p','Search','Year'], ascending=[True,True,True,False])
fills = {'A': PatternFill('solid', fgColor='E2EFDA'), 'B': PatternFill('solid', fgColor='FFF2CC')}
for r, row in enumerate(log[cols].itertuples(index=False), 2):
    for c, v in enumerate(row, 1):
        x = sc.cell(r, c, v); x.font = N; x.alignment = Alignment(vertical='top', wrap_text=(c in (5, 11)))
        if c == 7 and v: x.hyperlink = f'https://doi.org/{v}'; x.font = Font(name=F, size=10, color='0563C1', underline='single')
    p = row[9]
    if p in fills:
        for c in range(1, len(cols)+1):
            sc.cell(r, c).fill = fills[p]
sc.auto_filter.ref = f'A1:{get_column_letter(len(cols))}{len(log)+1}'
last = len(log) + 1
# Search log
sl = wb.create_sheet('Search log', 1)
header(sl, ['Search','Database','Date run','Exact search string','Hits','Kept (live)','Kept A (live)','Purpose'], [8,10,11,90,7,10,11,40])
purpose = {'S1':'Prior work: TF activity + chromatin in TNBC','S2':'Novelty of H3: TNBC TFs as CRISPR dependencies','S3':'TF activity as drug-response biomarker','S4':'Pharmacogenomic dataset consistency','S5':'ATAC-seq reproducibility (H1)','S6':'Regulon-method benchmarks (H2)','S7':'mTOR inhibitors in TNBC','S8a':'C-index with small samples / permutation','S8b':'PDX sample size / design','S9':'Cancer-biology reproducibility','S10':'TNBC cell-line subtypes'}
for i, s in enumerate(['S1','S2','S3','S4','S5','S6','S7','S8a','S8b','S9','S10'], 2):
    vals = [s, 'PubMed', '2026-11-05', strings.get(s, ''), counts[s]]
    for c, v in enumerate(vals, 1):
        x = sl.cell(i, c, v); x.font = N; x.alignment = Alignment(wrap_text=(c == 4), vertical='top')
    sl.cell(i, 6, f'=SUMPRODUCT(ISNUMBER(SEARCH("|"&A{i}&"|","|"&Screening!$C$2:$C${last}&"|"))*(Screening!$I$2:$I${last}="keep")*(Screening!$B$2:$B${last}="PubMed"))').font = N
    sl.cell(i, 7, f'=SUMPRODUCT(ISNUMBER(SEARCH("|"&A{i}&"|","|"&Screening!$C$2:$C${last}&"|"))*(Screening!$I$2:$I${last}="keep")*(Screening!$J$2:$J${last}="A")*(Screening!$B$2:$B${last}="PubMed"))').font = N
    sl.cell(i, 8, purpose[s]).font = N
tr = 13
sl.cell(tr, 1, 'Total').font = B
sl.cell(tr, 5, '=SUM(E2:E12)').font = B
sl.cell(tr+1, 1, 'Unique PubMed records').font = B
sl.cell(tr+1, 5, f'=COUNTIFS(Screening!$B$2:$B${last},"PubMed")').font = B
sl.cell(tr+2, 1, 'Kept overall (all sources)').font = B
sl.cell(tr+2, 5, f'=COUNTIFS(Screening!$I$2:$I${last},"keep")').font = B
sl.cell(tr+3, 1, 'Note: S8a needs a "C-index" re-run; its 98 hits gave only 1 useful record (methods papers were added by hand instead).  A record found by several searches is counted under each.').font = Font(name=F, size=9, italic=True)
# Key papers
kp = wb.create_sheet('Key papers')
kcols = ['Priority','Paper section','Citation','Title','DOI','Reason / role','Source']
header(kp, kcols, [8,30,32,70,30,55,14])
keep = log[log.Decision == 'keep'].sort_values(['_p','Paper section'])
for r, row in enumerate(keep[kcols].itertuples(index=False), 2):
    for c, v in enumerate(row, 1):
        x = kp.cell(r, c, v); x.font = N; x.alignment = Alignment(vertical='top', wrap_text=(c in (4, 6)))
        if c == 5 and v: x.hyperlink = f'https://doi.org/{v}'; x.font = Font(name=F, size=10, color='0563C1', underline='single')
    if row[0] in fills:
        for c in range(1, len(kcols)+1): kp.cell(r, c).fill = fills[row[0]]
kp.auto_filter.ref = f'A1:G{len(keep)+1}'
wb.save(f'{lit}/literature_log.xlsx'); print('saved', len(log), 'rows;', len(keep), 'kept')
