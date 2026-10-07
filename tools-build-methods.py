# Builds web/methods.html from web/MATH.md so the two never drift.
import re, html
src=open('web/MATH.md').read().split('\n')
out=[]; i=0; inlist=False; para=[]
def inline(t):
    t=html.escape(t,quote=False)
    t=re.sub(r'(https?://[^\s<]+)',r'<a href="\1">\1</a>',t)
    return t
def flush():
    global para
    if para: out.append('<p>'+inline(' '.join(para))+'</p>'); para=[]
while i<len(src):
    l=src[i]
    if l.startswith('\\['):
        flush(); blk=[]
        while not src[i].startswith('\\]'): blk.append(src[i]); i+=1
        blk.append(src[i]); out.append('<div class="math">'+html.escape('\n'.join(blk),quote=False)+'</div>')
    elif l.startswith('# '): flush(); out.append('<h1>'+inline(l[2:])+'</h1>')
    elif l.startswith('## '):
        flush()
        if inlist: out.append('</ul>'); inlist=False
        out.append('<h2>'+inline(l[3:])+'</h2>')
    elif re.match(r'^(- |\d+\. )',l):
        flush()
        if not inlist: out.append('<ul>'); inlist=True
        out.append('<li>'+inline(re.sub(r'^(- |\d+\. )','',l))+'</li>')
    elif l.strip()=='':
        flush()
        if inlist: out.append('</ul>'); inlist=False
    else: para.append(l)
    i+=1
flush()
body='\n'.join(out).replace('<h1>Censor Check: mathematical contract</h1>','<h1>The censored record,<br>with its assumptions showing.</h1>',1)
page=f'''<!doctype html><html lang="en" data-theme="dark"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Censor Check · Math, assumptions &amp; evidence · v0.1.0-alpha</title><link href="https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Instrument+Serif&display=swap" rel="stylesheet"><link rel="stylesheet" href="style.css">
<script>window.MathJax={{tex:{{inlineMath:[['\\\\(','\\\\)']],displayMath:[['\\\\[','\\\\]']]}},svg:{{fontCache:'global'}}}};</script><script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script></head><body>
<header><a href="index.html">← Censor Check</a><nav aria-label="Resources"><a href="model.js" download>Download model code</a><a href="MATH.md">Math as Markdown</a></nav></header>
<main class="methods"><span class="badge">DATA-RICH, INSIGHT-POOR · THE AUDIT, NO. 1 · 99 SMALL PROBLEMS · NO. 06 · v0.1.0-alpha</span>
{body}
</main></body></html>'''
open('web/methods.html','w').write(page); print('ok')
