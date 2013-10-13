def comb(axes, JPos, modulo, y=0.9, height=0.03, overhang=True):
    bottom, top = axes.get_ylim()
    left, right = axes.get_xlim()
    xlen = right - left
    ydiv = top - bottom
    textoffset = 0
    h = (top-bottom)*(y+height) + bottom
    for j, Pos in JPos:
        axes.axvline(x=Pos, ymin=y-height, ymax=y+height, color='k')
        if not j%modulo:
            if overhang:
                axes.axvline(x=Pos, ymin=y-height, ymax=y+height+0.04,
                             color='k')
                textoffset = 0.0025*xlen
            axes.text(Pos+textoffset, h+0.01*ydiv, r'%d' %j)
    axes.axhline(h, (JPos[0][1]-left)/xlen, (JPos[-1][1]-left)/xlen, color='k')
