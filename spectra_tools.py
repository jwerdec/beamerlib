def comb(axes, JPos, modulo, y=0.9, height=0.03, overhang=True, linewidth=0.5,
         extend=None, shiftnum=0):
    bottom, top = axes.get_ylim()
    left, right = axes.get_xlim()
    xlen = right - left
    ydiv = top - bottom
    hline_left = JPos[0][1]
    hline_right = JPos[-1][1]
    if extend == 'left' or extend == 'both':
        hline_left = left
    if extend == 'right' or extend == 'both':
        hline_right = right
    textoffset = 0
    h = (top-bottom)*(y+height) + bottom
    for j, Pos in JPos:
        axes.axvline(x=Pos, ymin=y-height, ymax=y+height, color='k',
                     linewidth=linewidth)
        if not (j+shiftnum)%modulo:
            if overhang:
                axes.axvline(x=Pos, ymin=y-height, ymax=y+height+0.04,
                             color='k', linewidth=linewidth)
                textoffset = 0.0025*xlen
            axes.text(Pos+textoffset, h+0.01*ydiv, r'%d' %j)
    axes.axhline(h, (hline_left-left)/xlen, (hline_right-left)/xlen, color='k',
                 linewidth=linewidth)
