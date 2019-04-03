#!/usr/local/sci/bin/python2.7
'''
@author: Andrew Ciavarella
Mod. of qfigs_0.py to test display of p0,p1 bounds.
Produce quick figures for event attribution.
'''

import matplotlib.pyplot as plt
plt.rcdefaults()

import matplotlib.patches as mpatches
import numpy as np

def rects(xy, w, p0,p1,cols,alpha=None,ec=None,lw=None,overlap=False):
    # ,p0_min=None,p0_max=None,p1_min=None,p1_max=None
    '''Return a pair (ALL & NAT) of rectangle patches'''

    # NAT is shifted right by a small amount from ALL
    if overlap:
        xy1 = [xy[0]+0.45*w, xy[1]]
        xy0 = [xy[0]+0.55*w, xy[1]]
    else:
        xy1 = xy
        xy0 = [xy[0]+1.0*w, xy[1]]
    rect_p1 = mpatches.Rectangle(xy1, w, p1, color=cols['all'], ec=ec, alpha=alpha, lw=lw)
    rect_p0 = mpatches.Rectangle(xy0, w, p0, color=cols['nat'], ec=ec, alpha=alpha, lw=lw)

    return rect_p1, rect_p0    

def barchart_p1p0(events, bounds=False):
    '''
    Return figure and axis objects containing simple bar chart of multiple event p1, p0 values.
    Argument 'events' is list of Event instances, defined in opatt/event_def.py.
    Bar chart should always look 'nice': code written so that labels are respectful of bars and edges of plot.
    '''

    # Color palette: CMIP5 D&A colours from Fraser Lott.
    cols = {'all':'#c80000', 'nat':'#00b478', 'ant':'#b8890a', 'ghg':'#4646ff', 'oan':'#ffc800'}

    # Width of region column etc: depends on no. regions
    n = len(events)
    dx = 1.0/float(n)
    # Width of bars
    w = 0.4*dx
    # Positions of lower left of bars
    xs = [(0.1+float(i))*dx for i in range(n)]
    ys = [0.0 for i in range(n)]

    # Create figure
    fig, ax = plt.subplots()

    # X limit leave room either side. Fixed interval for any no. regions.
    plt.xlim(-0.1,1.1)


    # Y limit leave room for legend text
    # Adjust max y to account for small max p1, p0 values
    p_list = []
    for event in events:
        if not bounds:
            try:
                p_list.append(event.p0[1])
                p_list.append(event.p1[1])
            except TypeError:#means we only have the central probability
                p_list.append(event.p0)
                p_list.append(event.p1)
            except IndexError:
                p_list.append(event.p0)
                p_list.append(event.p1)
        else:
            p_list.extend(event.p0)
            p_list.extend(event.p1)

    p_max = max(p_list)
    pc_max = 100.0*p_max
    y_min = -15.0
    if pc_max > 10.0:
        y_max = 10.0*(int(pc_max/10.0)+3)
    else:
        if pc_max > 1.0:
            y_max = int(pc_max)+3
            y_min = -0.3333*y_max
        else:
            y_max = 0.1*(int(pc_max*10.0)+3)
            y_min = -0.5*y_max

    plt.ylim(y_min,y_max)
    #print '--- y_max = ', y_max
    #print '   pc_max = ', pc_max
    #print '    p_max = ', p_max

    # Generate list of bar-pair patches (list of 2-element lists)
    patches = []
    for i, event in enumerate(events):
        # Take best estimate event probs
        try:
            p0 = event.p0[1]
            p1 = event.p1[1]
        except TypeError:#p1 and p0 and a single value, i.e. no error bars
            p0 = event.p0
            p1 = event.p1
        except IndexError:
            p0 = event.p0
            p1 = event.p1
        pc0 = 100.0*p0
        pc1 = 100.0*p1
        xy = [xs[i], ys[i]]
        alpha = 1.0
        ec = None
        lw = None
        if bounds:
            alpha = 0.5
            ec = "black"
            lw = 2.0
        rect_p1, rect_p0 = rects(xy,w,pc0,pc1,cols,alpha=alpha,ec=ec,lw=lw,
                                 overlap=not bounds)

        if not bounds:
            rects_list = [rect_p1, rect_p0]
        else:
            fdx = 0.0
            xy_min = [xs[i]+fdx*w, ys[i]]
            rect_p1_min, rect_p0_min = rects(xy_min,w,100.0*event.p0[0],
                                            100.0*event.p1[0],cols,alpha=1.0,ec="none")
            xy_max = [xs[i]-fdx*w, ys[i]]
            rect_p1_max, rect_p0_max = rects(xy_max,w,100.0*event.p0[2],
                                             100.0*event.p1[2],cols,alpha=0.2,ec="none")            
            rects_list = [rect_p1, rect_p0,
                          rect_p1_min, rect_p0_min,
                          rect_p1_max, rect_p0_max]

        patches.append(rects_list)

    # Add bars to figure
    for i in range(n):
        if not bounds:
            ax.add_patch(patches[i][0])
            ax.add_patch(patches[i][1])
        else:
            ax.add_patch(patches[i][2])
            ax.add_patch(patches[i][0])
            ax.add_patch(patches[i][4])
            ax.add_patch(patches[i][3])
            ax.add_patch(patches[i][1])
            ax.add_patch(patches[i][5])            

    # Region labels below bars
    for i, event in enumerate(events):
        plt.text(xs[i]+0.5*dx, 0.66*y_min, event.reg, ha="right", family='sans-serif', size=12)

    # Y-axis label
    plt.ylabel('% Chance of occurrence')

    # Truncate Y-axis ticks
    yticks = ax.get_yticks()
    yticks_trunc = [yt for yt in yticks if ((yt <= 100.0) and (yt >= 0.0))]
    ax.set_yticks(yticks_trunc)

    # Remove all X-axis ticks
    ax.set_xticks([])

    # Legend text: place at height respectful of bars and top-most gridline
    if y_max > 100.0:
        if y_max > 110.0:
            txty1 = 112.0
            txty0 = 103.0
        else:
            yran = y_max - pc_max
            txty1 = pc_max + 0.7*yran
            txty0 = pc_max + 0.3*yran
    else:
        dytick = yticks_trunc[-1] - yticks_trunc[-2]
        txty1 = yticks_trunc[-2] + 0.6*dytick
        txty0 = yticks_trunc[-2] + 0.2*dytick
    plt.text(0.1, txty1, 'With human influence', ha="left", family='sans-serif', size=14, color=cols['all'])
    plt.text(0.1, txty0, 'Without human influence', ha="left", family='sans-serif', size=14, color=cols['nat'])

    ax.yaxis.grid(True)
    ax.grid(which='major', axis='y', linestyle='dotted', color='gray', linewidth=2)
    ax.set_axisbelow(True)

    return fig, ax

    #plt.show()

####################### MAIN ###########################

def main(n):
    '''
    @author: Andrew Ciavarella
    Mod. of test_rects2.py.
    Produce quick figures for event attribution.
    '''

if __name__ == "__main__":
    main()
