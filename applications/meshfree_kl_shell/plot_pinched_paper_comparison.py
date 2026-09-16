#!/usr/bin/env python3
"""Create an SVG comparison from fresh pinched-cylinder CSV results.

The nonlinear paper curve is not digitized in this repository. Its documented
<1000 value at 150 mm is shown as an upper-bound marker, not an invented curve.
This script intentionally uses only the Python standard library.
"""
from pathlib import Path
import csv, math, struct, zlib

ROOT = Path(__file__).resolve().parent
with open(ROOT / "pinched_current_comparison.csv") as f:
    linear = list(csv.DictReader(f))
with open(ROOT / "pinched_plastic_current.csv") as f:
    plastic = list(csv.DictReader(f))

W, H = 1200, 500
panels = [(70, 55, 500, 370), (660, 55, 500, 370)]
parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
         '<rect width="100%" height="100%" fill="white"/>',
         '<style>text{font-family:sans-serif;fill:#222}.axis{stroke:#222}.grid{stroke:#ddd}.paper{stroke:#555;stroke-dasharray:7 5}.old{stroke:#c44}.new{stroke:#1764aa}.coll{stroke:#e28b16}</style>']

def axes(box, title, xlabel, ylabel, xmin, xmax, ymin, ymax):
    x,y,w,h=box
    for i in range(6):
        xx=x+i*w/5; yy=y+i*h/5
        parts.append(f'<line class="grid" x1="{xx}" y1="{y}" x2="{xx}" y2="{y+h}"/>')
        parts.append(f'<line class="grid" x1="{x}" y1="{yy}" x2="{x+w}" y2="{yy}"/>')
        parts.append(f'<text x="{xx}" y="{y+h+20}" font-size="11" text-anchor="middle">{xmin+(xmax-xmin)*i/5:.0f}</text>')
        parts.append(f'<text x="{x-8}" y="{y+h-i*h/5+4}" font-size="11" text-anchor="end">{ymin+(ymax-ymin)*i/5:.1f}</text>')
    parts.append(f'<rect class="axis" fill="none" x="{x}" y="{y}" width="{w}" height="{h}"/>')
    parts.append(f'<text x="{x+w/2}" y="{y-22}" font-size="17" text-anchor="middle">{title}</text>')
    parts.append(f'<text x="{x+w/2}" y="{y+h+45}" font-size="13" text-anchor="middle">{xlabel}</text>')
    parts.append(f'<text transform="translate({x-52},{y+h/2}) rotate(-90)" font-size="13" text-anchor="middle">{ylabel}</text>')
    return lambda X,Y:(x+(X-xmin)/(xmax-xmin)*w, y+h-(Y-ymin)/(ymax-ymin)*h)

def poly(points, mapper, cls):
    p=' '.join(f'{mapper(a,b)[0]:.1f},{mapper(a,b)[1]:.1f}' for a,b in points)
    parts.append(f'<polyline class="{cls}" fill="none" stroke-width="2" points="{p}"/>')
    for a,b in points:
        x,y=mapper(a,b); parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3" fill="currentColor" class="{cls}"/>')

m=axes(panels[0], 'Linear pinched cylinder', 'nodes', 'displacement / paper', 0, 9000, 0, 7)
for mode,cls in [('deck_before_fix_cubic','old'),('paper_quadratic','new')]:
    d=[(float(r['nodes']),float(r['ratio_to_paper'])) for r in linear if r['configuration']==mode]
    poly(d,m,cls)
x1,y1=m(0,1);x2,y2=m(9000,1);parts.append(f'<line class="paper" x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}"/>')
parts += ['<text x="90" y="82" class="old" font-size="12">old cubic-completeness deck</text>',
          '<text x="90" y="100" class="new" font-size="12">paper quadratic configuration</text>',
          '<text x="90" y="118" font-size="12">-- paper reference 1.8248e-5</text>']

m=axes(panels[1], 'Elasto-plastic pinch (nt=40)', 'physical displacement (mm)', 'raw generalized reaction', 0, 160, 0, 35000)
for mode,cls in [('coefficient_KBC','new'),('collocation_KBC','coll')]:
    d=[(float(r['physical_displacement_mm']),float(r['reported_generalized_reaction'])) for r in plastic if r['bc_mode']==mode]
    poly(d,m,cls)
px,py=m(150,1000);parts.append(f'<path d="M {px-6} {py-8} L {px+6} {py-8} L {px} {py+4} Z" fill="#222"/>')
parts += ['<text x="680" y="82" class="new" font-size="12">coefficient KBC (raw reaction)</text>',
          '<text x="680" y="100" class="coll" font-size="12">collocation KBC (raw reaction)</text>',
          f'<text x="{px-8}" y="{py-15}" font-size="11" text-anchor="end">paper: &lt;1000 at 150 mm (upper bound)</text>',
          '<text x="680" y="410" font-size="11">Collocation generalized reactions require a physical-force transform.</text>']
parts.append('</svg>')
(ROOT/'pinched_paper_comparison.svg').write_text('\n'.join(parts))
print(ROOT/'pinched_paper_comparison.svg')

# Also emit a dependency-free raster version for reports and review systems that
# do not render SVG. This deliberately small rasterizer draws the same data,
# axes, tick marks, and a compact bitmap-font legend.
PW, PH = 1200, 500
pixels = bytearray([255, 255, 255] * PW * PH)
COLORS = {'black':(30,30,30), 'grid':(220,220,220), 'old':(196,68,68),
          'new':(23,100,170), 'coll':(226,139,22), 'paper':(85,85,85)}

def pixel(x, y, color):
    x, y = int(x), int(y)
    if 0 <= x < PW and 0 <= y < PH:
        k = 3*(y*PW+x); pixels[k:k+3] = bytes(COLORS[color])

def line(x0, y0, x1, y1, color, width=1, dashed=False):
    x0,y0,x1,y1=map(int,(x0,y0,x1,y1)); dx=abs(x1-x0); sx=1 if x0<x1 else -1
    dy=-abs(y1-y0); sy=1 if y0<y1 else -1; err=dx+dy; n=0
    while True:
        if not dashed or (n//7)%2 == 0:
            for ox in range(-(width//2), width//2+1):
                for oy in range(-(width//2), width//2+1): pixel(x0+ox,y0+oy,color)
        if x0==x1 and y0==y1: break
        e2=2*err
        if e2>=dy: err+=dy; x0+=sx
        if e2<=dx: err+=dx; y0+=sy
        n+=1

def dot(x,y,color,r=3):
    for yy in range(-r,r+1):
        for xx in range(-r,r+1):
            if xx*xx+yy*yy<=r*r: pixel(x+xx,y+yy,color)

FONT = {
' ':['000','000','000','000','000'],'-':['000','000','111','000','000'],'.':['000','000','000','000','010'],
'0':['111','101','101','101','111'],'1':['010','110','010','010','111'],'2':['111','001','111','100','111'],
'3':['111','001','111','001','111'],'4':['101','101','111','001','001'],'5':['111','100','111','001','111'],
'6':['111','100','111','101','111'],'7':['111','001','010','010','010'],'8':['111','101','111','101','111'],
'9':['111','101','111','001','111'],'A':['010','101','111','101','101'],'B':['110','101','110','101','110'],
'C':['111','100','100','100','111'],'D':['110','101','101','101','110'],'E':['111','100','110','100','111'],
'F':['111','100','110','100','100'],'G':['111','100','101','101','111'],'H':['101','101','111','101','101'],
'I':['111','010','010','010','111'],'J':['001','001','001','101','111'],'K':['101','101','110','101','101'],
'L':['100','100','100','100','111'],'M':['101','111','111','101','101'],'N':['101','111','111','111','101'],
'O':['111','101','101','101','111'],'P':['111','101','111','100','100'],'Q':['111','101','101','111','001'],
'R':['110','101','110','101','101'],'S':['111','100','111','001','111'],'T':['111','010','010','010','010'],
'U':['101','101','101','101','111'],'V':['101','101','101','101','010'],'W':['101','101','111','111','101'],
'X':['101','101','010','101','101'],'Y':['101','101','010','010','010'],'Z':['111','001','010','100','111'],
':':['000','010','000','010','000'],'/':['001','001','010','100','100'],'<':['001','010','100','010','001'],
'=':['000','111','000','111','000'],'%':['101','001','010','100','101'],'_':['000','000','000','000','111']}

def text(x,y,s,color='black',scale=2):
    ox=x
    for ch in s.upper():
        glyph=FONT.get(ch,FONT[' '])
        for gy,row in enumerate(glyph):
            for gx,v in enumerate(row):
                if v=='1':
                    for a in range(scale):
                        for b in range(scale): pixel(x+gx*scale+a,y+gy*scale+b,color)
        x += 4*scale
    return x-ox

def raster_axes(box, xmin,xmax,ymin,ymax):
    x,y,w,h=box
    for i in range(6):
        xx=x+i*w/5; yy=y+i*h/5
        line(xx,y,xx,y+h,'grid'); line(x,yy,x+w,yy,'grid')
        text(xx-10,y+h+8,f'{xmin+(xmax-xmin)*i/5:.0f}',scale=1)
        text(x-34,y+h-i*h/5-2,f'{ymin+(ymax-ymin)*i/5:.1f}',scale=1)
    line(x,y,x+w,y,'black');line(x,y+h,x+w,y+h,'black');line(x,y,x,y+h,'black');line(x+w,y,x+w,y+h,'black')
    return lambda X,Y:(x+(X-xmin)/(xmax-xmin)*w, y+h-(Y-ymin)/(ymax-ymin)*h)

def raster_poly(data, mapper, color):
    pts=[mapper(a,b) for a,b in data]
    for a,b in zip(pts,pts[1:]): line(*a,*b,color,width=2)
    for p in pts: dot(*p,color)

rm=raster_axes(panels[0],0,9000,0,7)
text(170,18,'LINEAR PINCHED CYLINDER',scale=2)
for mode,color in [('deck_before_fix_cubic','old'),('paper_quadratic','new')]:
    raster_poly([(float(r['nodes']),float(r['ratio_to_paper'])) for r in linear if r['configuration']==mode],rm,color)
a=rm(0,1);b=rm(9000,1);line(*a,*b,'paper',dashed=True)
line(90,72,120,72,'old',2);text(130,67,'OLD CUBIC DECK','old',1)
line(90,88,120,88,'new',2);text(130,83,'PAPER QUADRATIC','new',1)
line(90,104,120,104,'paper',1,True);text(130,99,'PAPER REFERENCE','paper',1)
text(270,465,'NODES',scale=2);text(76,438,'RATIO TO PAPER',scale=1)

rm=raster_axes(panels[1],0,160,0,35000)
text(770,18,'ELASTO PLASTIC PINCH NT 40',scale=2)
for mode,color in [('coefficient_KBC','new'),('collocation_KBC','coll')]:
    raster_poly([(float(r['physical_displacement_mm']),float(r['reported_generalized_reaction']))
                 for r in plastic if r['bc_mode']==mode],rm,color)
x,y=rm(150,1000);dot(x,y,'black',5)
line(680,72,710,72,'new',2);text(720,67,'COEFFICIENT KBC RAW','new',1)
line(680,88,710,88,'coll',2);text(720,83,'COLLOCATION KBC RAW','coll',1)
dot(695,104,'black',4);text(720,99,'PAPER < 1000 AT 150 MM','black',1)
text(815,465,'PHYSICAL DISPLACEMENT MM',scale=2)

def chunk(tag, data):
    return struct.pack('>I',len(data))+tag+data+struct.pack('>I',zlib.crc32(tag+data)&0xffffffff)
raw=b''.join(b'\x00'+pixels[y*PW*3:(y+1)*PW*3] for y in range(PH))
png=b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',PW,PH,8,2,0,0,0))+chunk(b'IDAT',zlib.compress(raw,9))+chunk(b'IEND',b'')
(ROOT/'pinched_paper_comparison.png').write_bytes(png)
print(ROOT/'pinched_paper_comparison.png')
