#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Experimental functions for IPE beamline"""


# %% =========================== fancy plot =============================== %% #
def _plotall(filepath, scan, axes=None, set_window_size=True):
    TEY, TFY, I0, PD = read(filepath, scan) 

    if axes is None:
        fig, axes = br.subplots(2, 3)
    if set_window_size:
        br.set_window_size((514, 1091)) 
        plt.subplots_adjust(left=0.06, right=0.99, top=0.95, bottom=0.1)

    seq = [TEY, MCP, TFY, TEY/RMU, MCP/RMU, TFY/RMU]
    seq[3]. label = 'TEY/RMU'
    seq[4]. label = 'MCP/RMU'
    seq[5]. label = 'TFY/RMU'
    # seq[6]. label = 'norm TEY/RMU'
    # seq[7]. label = 'norm MCP/RMU'
    # seq[8]. label = 'norm TFY'
    for i in range(2*3):
        ax = axes[i]
        seq[i].plot(ax=ax, label=f'{seq[i].label}, #{seq[i].scan}, x={round(seq[i].sample_x, 4)}, y={round(seq[i].sample_y, 4)}')
        ax.labels_xas()
        br.leg(ax=ax, fontsize='xx-small')
    return axes

def _sequential(*args, **kwargs):

    ################
    # process data #
    ################
    temp = _process(folderpath=folderpath, sbins=sbins, calib=calib, norm=norm)
    s    = temp['s']
    pe1  = temp['pe1']
    pe2  = temp['pe2']
    pes1 = temp['pes1']
    pes2 = temp['pes2']

    #######################
    # initial definitions #
    #######################
    pes1.__i  = 0

    ######################
    # change keybindings #
    ######################
    try:
        matplotlib.rcParams['keymap.back'].remove('left')
        matplotlib.rcParams['keymap.forward'].remove('right')
    except ValueError:
        pass

    ###################
    # keyboard events #
    ###################
    def keyboard(event, pes1, pes2, axes):
        if event.key == 'right':
            # increase i
            pes1.__i = pes1.__i + 1
            if pes1.__i >= len(pes1):
                pes1.__i = len(pes1) - 1
    
        elif event.key == 'left':# or event.key == 'down':
            # decrease i
            pes1.__i = pes1.__i - 1
            if pes1.__i < 0:
                pes1.__i = 0
        else:
            return
            
        # clear axis
        axes[0].cla()
        axes[1].cla()
        
        # set labels
        axes[0].set_xlabel('x (pixel)')
        axes[0].set_ylabel('y (pixel)')
        axes[1].set_xlabel('counts/bin')
        
        # change title
        axes[0].set_title('Use left/right keyboard keys to flip through images: ' + str(pes1.__i) + '/' + str(len(pes1)-1), fontsize='small')

        # plot axes 0
        pes1[pes1.__i].plot(ax=axes[0], show_limits=True, **kwargs)
        pes2[pes1.__i].plot(ax=axes[0], show_limits=True, **kwargs)

        # plot axes 1
        pes1[pes1.__i].integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[1])
        pes2[pes1.__i].integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[1])
    
        plt.draw()

    #######################
    # initializing figure #
    #######################
    fig, axes = br.subplots(4, 2, width_ratios=[4, 1], height_ratios=[1, 1, 1, 2], wspace=0.1, hspace=0.8, figsize=(18, 26))
    axes[1].remove_yticklabels()
    axes[3].remove_yticklabels()
    
    ##############
    # share axis #
    ##############
    br.sharey([axes[0], axes[1]])
    br.sharey([axes[2], axes[3]])
    

    ##################
    # error messages #
    ##################
    if pe1.RIXSCam_NumImages != len(pes1):
        fig.suptitle(f'ERROR: # of images ({len(pes1)}) inside folder is different from # of acquired images ({int(pe1.RIXSCam_NumImages)})', color='red')

    ######################
    # set initial titles #
    ######################
    axes[0].set_title('Use left/right keyboard keys to flip through images: ' + str(0) + '/' + str(len(pes1)-1), fontsize='small')
    axes[1].set_title(f'nbins = {sbins}', fontsize='small')
    axes[2].set_title('Summed photon events for each CCD', fontsize='small')
    axes[4].set_title('Number of photons per image', fontsize='small')
    axes[6].set_title('Final spectrum', fontsize='small')
    

    ########
    # plot #
    ########
    # plot initial photon events (axes 0)
    pes1[0].plot(ax=axes[0], show_limits=True, **kwargs)
    pes2[0].plot(ax=axes[0], show_limits=True, **kwargs)

    # plot initial spectra (axes 1)
    pes1[0].integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[1])
    pes2[0].integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[1])

    # plot photon events summed (axes 2)
    pe1.plot(ax=axes[2], show_limits=True, **kwargs)
    pe2.plot(ax=axes[2], show_limits=True, **kwargs)

    # plot spectra summed (axes 3)
    pe1.integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[3])
    pe2.integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[3])

    # plot number of photons per image (axes 4)
    for pes in (pes1, pes2):
        number_of_photons_ccd = [len(_pe) for _pe in pes]
        axes[4].plot(np.arange(0, len(number_of_photons_ccd)), number_of_photons_ccd, marker='o', lw=1)

    # plot spectrum (axes 6)
    s.plot(ax=axes[6], color='black')

    ##############
    # set labels #
    ##############
    for i in (0, 2):
        axes[i].set_xlabel('x (pixel)')
        axes[i].set_ylabel('y (pixel)')
    for i in (1, 3):
        axes[i].set_xlabel('counts/bin')
    axes[4].set_xlabel('Image number')
    axes[4].set_ylabel('Number of photons')

    if calib is None:
        axes[6].set_xlabel('y (pixel)')
    else:
        axes[6].set_xlabel('Energy (eV)')
    if norm:
        axes[6].set_ylabel('Norm. intensity (arb. units)')
    else:
        axes[6].set_ylabel('Photon count per bin')

    ######################
    # register callbacks #
    ######################
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, pes1=pes1, pes2=pes2, axes=axes))
    return temp

# from brixs.addons.png2clipboard import png2clipboard
# def figure2clipboard():
#     plt.savefig(TMP/'temp.png', dpi=1000)
#     png2clipboard(TMP/'temp.png')
