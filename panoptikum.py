# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 16:12:46 2013

@author: Reinhardt A.W. Maier <rema@zaehlwerk.net>
"""
import gtk
import apogee
import loadColormap
import imageAnalyse
import math
import gio
import os
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as Nav
from matplotlib.widgets import RectangleSelector
from numpy import size, array, linspace, sort, zeros, abs, shape
from pylab import random

class Panoptikum(object):
    def __init__(self):
        # assign gtk-Layout
        builder = gtk.Builder()
        builder.add_from_file('panoptikum.glade')
        builder.connect_signals(self)
        self.window = builder.get_object('window1')
        self.windowMessage = builder.get_object('messagedialog1')
        self.filechooserProject = builder.get_object('filechooserdialog1')
        self.filefilter = builder.get_object('filefilter1')
        self.filefilter.add_pattern("*.cfg")
        self.large_atom_check_button = builder.get_object('large_atom_check_button')
        self.live_mode_check_button = builder.get_object('live_mode_check_button')
        self.roiOnly_check_button = builder.get_object('roiOnly_check_button')
        self.sameID_check_button = builder.get_object('sameID_check_button')
        self.radio_buttons_imageCategory = builder.get_object('action1')
        self.window.connect('destroy', self.__del__)
        self.statusbar = builder.get_object('statusbar1')
        self.hbox_Rb = builder.get_object('hboxRb')
        self.vbox_Rb = builder.get_object('vboxRb')
        self.vbox_RbLast10 = builder.get_object('vboxRbLast10')
        self.hbox_Li = builder.get_object('hboxLi')
        self.vbox_Li = builder.get_object('vboxLi')
        self.www_Rb = builder.get_object('togglebuttonRbWww')
        self.www_Li = builder.get_object('togglebuttonLiWww')
        
        # global Variables
        self.MessageType = 'N_Li'
        self.maxImageSize = array([0, 1392, 0, 1040])

        ######### TODO ############
        """
        self.TEST = builder.get_object('filechooserbutton1Rb')
        self.TEST2 = builder.get_object('entry1')
        #import time
        #self.time_start = 0
        #self.time_end = 0
        try:
            print p['pixelSizeTextField'].get_value()
            print self.TEST.get_filename()
            self.TEST.set_filename('rb_recenta.fit')
            self.TEST2.set_text(str(p['part'])[1:-1])
        except:
            pass
        """
        """
        self.time_start = time.clock()
                self.time_end = time.clock()
        print self.time_end - self.time_start
        """
        ######### TODO ############

        # matplotlib
        self.whiteCMap = loadColormap.whiteCMap()
        self.initNav()
        self.shiftPressed = False

        # Rb __init__
        self.figure_Rb = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
        ax1 = plt.subplot(gs[0])
        ax1.set_xlabel(u'\u00b5m')
        ax1.set_ylabel(u'\u00b5m')
        ax1.xaxis.set_label_coords(0.97, -0.06)
        ax3 = plt.subplot(gs[1])
        ax3.xaxis.set_label_coords(0.97, -0.16)
        canvas_Rb = FigCanvas(self.figure_Rb)
        plt.subplots_adjust(hspace=0.25, bottom=0.05, right=0.95 )
        nav_Rb = Nav(canvas_Rb, self.window)
        self.vbox_Rb.pack_start(nav_Rb, False, False)
        self.vbox_Rb.pack_start(canvas_Rb)
        self.figure_RbLast10 = plt.figure()
        self.ax1Last10 = plt.subplot(111)
        self.canvas_RbLast10 = FigCanvas(self.figure_RbLast10)
        self.vbox_RbLast10.pack_start(self.canvas_RbLast10)
        RbValues2calc = ['N_Rb','T_Rb','fwhmV_Rb','fwhmH_Rb','PosV_Rb','PosH_Rb','AmpV_Rb','AmpH_Rb','OffV_Rb','OffH_Rb','______AUTOMATISIERUNG_','AUTO_N_POS','AUTO_STEP_POS','AUTO_START','AUTO_STEP','AUTO_BREAK','ABB_TOF_RB']
        RbIndexLst = [s.lower().strip('_' + 'rb') for s in RbValues2calc]
        self.RbParameters = {
            'pixelSize' : 6.57e-6,
            'pixelSizeTextField' : builder.get_object('spinbutton1'),
            'axScale' : 1e6,
            'species' : 'Rb',
            'folder' : 'F:\\INBOX\\Apogee\\',
            'filenames' : ['rb_recenta.fit', 'rb_recentb.fit', 'rb_recentc.fit'],
            'adwinfile' : '..\\ADwin\\adwin_code_recent.acm',
            'savefilename' : 'rb_result.png',
            'fullFrame' : zeros(2),
            'partAtomCount' : array([550, 750, 500, 800]),
            'partLumCorr' : array([550, 750, 800, 950]),
            'part' : self.maxImageSize,
            'xylim' : [(0, 6826), (9138, 0)],
            'Marker1Pos': (0,0),
            'Marker2Pos': (0,0),
            'imageAx' : ax1,
            'linescanAx' : ax3,
            'canvasObj' : canvas_Rb,
            'ODmaxAuto' : builder.get_object('checkbuttonRbAutoOD'),
            'ODmaxLabel' : builder.get_object('labelRbOD'),
            'ODmax' : builder.get_object('hscaleRbOD'),
            'adwinID' : builder.get_object('adwinIDRb'),
            'imageCategory': 'opt. Dichte',
            'values2calc' : RbValues2calc,
            'indexLst': RbIndexLst,
            'last10' : pd.DataFrame(index = RbIndexLst)
            }

        # Li __init__
        self.figure_Li = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
        ax0 = plt.subplot(gs[0])
        ax0.set_xlabel(u'\u00b5m')
        ax0.set_ylabel(u'\u00b5m')
        ax0.xaxis.set_label_coords(0.97, -0.06)
        ax2 = plt.subplot(gs[1])
        ax2.xaxis.set_label_coords(0.97, -0.16)
        canvas_Li = FigCanvas(self.figure_Li)
        plt.subplots_adjust(hspace=0.25, bottom=0.05, right=0.95 )
        nav_Li = Nav(canvas_Li, self.window)
        self.vbox_Li.pack_start(nav_Li, False, False)
        self.vbox_Li.pack_start(canvas_Li)
        LiValues2calc = ['N_Li','T_Li','fwhmV_Li','fwhmH_Li','PosV_Li','PosH_Li','AmpV_Li','AmpH_Li','OffV_Li','OffH_Li','______AUTOMATISIERUNG_','AUTO_N_POS','AUTO_STEP_POS','AUTO_START','AUTO_STEP','AUTO_BREAK','ABB_TOF_LI']
        LiIndexLst = [s.lower().strip('_' + 'li') for s in LiValues2calc]
        self.LiParameters = {
            'pixelSize' : 6.57e-6,
            'axScale' : 1e6,
            'species' : 'Li',
            'folder' : 'F:\\INBOX\\Apogee\\',
            'filenames' : ['li_recenta.fit', 'li_recentb.fit', 'li_recentc.fit'],
            'adwinfile' : '..\\ADwin\\adwin_code_recent.acm',
            'savefilename' : 'li_result.png',
            'fullFrame' : zeros(2),
            'partAtomCount' : array([550, 750, 500, 800]),
            'partLumCorr' : array([550, 750, 800, 950]),
            'part' : self.maxImageSize,
            'xylim' : [(0, 6826), (9138, 0)],
            'Marker1Pos': (0,0),
            'Marker2Pos': (0,0),
            'imageAx' : ax0,
            'linescanAx' : ax2,
            'canvasObj' : canvas_Li,
            'ODmaxAuto' : builder.get_object('checkbuttonLiAutoOD'),
            'ODmaxLabel' : builder.get_object('labelLiOD'),
            'ODmax' : builder.get_object('hscaleLiOD'),
            'adwinID' : builder.get_object('adwinIDLi'),
            'imageCategory': 'opt. Dichte',
            'values2calc' : LiValues2calc,
            'indexLst': LiIndexLst,
            'last10' : pd.DataFrame(index = LiIndexLst)
            }

        # reload Parameters of last Session
        try:
            with open('restoreLastSession.pickle', 'rb') as infile:
                parameterFromStorage = pickle.load(infile)
                self.RbParameters['part'] = parameterFromStorage[0]
                self.RbParameters['partAtomCount'] = parameterFromStorage[1]
                self.RbParameters['partLumCorr'] = parameterFromStorage[2]
                self.RbParameters['xylim'] = parameterFromStorage[3]
                self.RbParameters['Marker1Pos'] = parameterFromStorage[4]
                self.RbParameters['Marker2Pos'] = parameterFromStorage[5]
                self.LiParameters['part'] = parameterFromStorage[6]
                self.LiParameters['partAtomCount'] = parameterFromStorage[7]
                self.LiParameters['partLumCorr'] = parameterFromStorage[8]
                self.LiParameters['xylim'] = parameterFromStorage[9]
                self.LiParameters['Marker1Pos'] = parameterFromStorage[10]
                self.LiParameters['Marker2Pos'] = parameterFromStorage[11]
        except:
            print 'Cannot load last Session!'

        # set working directory
        os.chdir(self.RbParameters['folder'])
        
        # draw matplotlib stuff
        self.initLines(self.RbParameters)
        self.initLines(self.LiParameters)
        self.updateFullFrameImage(self.RbParameters)
        self.updateFullFrameImage(self.LiParameters)
        self.last10Rb = self.drawImageArea(self.RbParameters)
        self.last10Li = self.drawImageArea(self.LiParameters)
        
        # show gtk-window
        self.window.set_default_size(1024, 800)
        self.window.set_size_request(600, 800)
        self.window.set_title ('Panoptikum - LiRb-Lab Image Analyse')
        self.window.show_all()
        #self.statusbar1_text_pushed('Li image updated... waiting for images...')

        # GTK-event handlers
        self.window.connect('key_press_event', self.on_key_press)
        self.window.connect('key_release_event', self.on_key_released)

        # matplotlib-event handlers
        canvas_Rb.mpl_connect('button_press_event', self.mouse_press_callback_Rb)
        canvas_Rb.mpl_connect('button_release_event', self.mouse_release_callback_Rb)
        canvas_Rb.mpl_connect('scroll_event', self.mouse_scrolled_callback_Rb)
        canvas_Li.mpl_connect('button_press_event', self.mouse_press_callback_Li)
        canvas_Li.mpl_connect('button_release_event', self.mouse_release_callback_Li)
        canvas_Li.mpl_connect('scroll_event', self.mouse_scrolled_callback_Li)
        rectprops = dict(facecolor = 'black', edgecolor = 'black', alpha = 0.1, fill = True)
        self.RS_Rb = RectangleSelector(ax1, self.RectangleSelector_callback_Rb,
                                       drawtype='box', useblit=True,
                                       button=[1,3], minspanx=10, minspany=10,
                                       spancoords='pixels', rectprops = rectprops)

        self.RS_Li = RectangleSelector(ax0, self.RectangleSelector_callback_Li,
                                       drawtype='box', useblit=True,
                                       button=[1,3], minspanx=10, minspany=10,
                                       spancoords='pixels', rectprops = rectprops)

        # GIO-event handlers
        fileRb = gio.File(self.RbParameters['filenames'][-1])
        self.monitorRb = fileRb.monitor_file()
        self.monitorRbID = self.monitorRb.connect('changed', self.file_changedRb)

        fileLi = gio.File('li_recentc.fit')
        self.monitorLi = fileLi.monitor_file()
        self.monitorLiID = self.monitorLi.connect('changed', self.file_changedLi)

    def initLines(self, parameters):
        p = parameters
        p['imageObj'] = p['imageAx'].imshow(random([2,2]), self.whiteCMap, interpolation='none', aspect='auto')
        p['lineFitH'], = p['imageAx'].plot(array([0, 1040])*p['pixelSize']*p['axScale'], [-1e5,-1e5], ':', lw=2, color = '0.15')
        p['lineFitV'], = p['imageAx'].plot([-1e5,-1e5], array([0, 1392])*p['pixelSize']*p['axScale'], ':', lw=2, color = '0.15')
        p['lineMarker1H'], = p['imageAx'].plot(array([0, 1040])*p['pixelSize']*p['axScale'], [-1e5,-1e5], 'k-', lw=1)
        p['lineMarker1V'], = p['imageAx'].plot([-1e5,-1e5], array([0, 1382])*p['pixelSize']*p['axScale'], 'k-', lw=1)
        self.drawCrosshair(p['lineMarker1H'], p['lineMarker1V'], p['Marker1Pos'])
        p['lineMarker2H'], = p['imageAx'].plot(array([0, 1040])*p['pixelSize']*p['axScale'], [-1e5,-1e5], 'r-', lw=1)
        p['lineMarker2V'], = p['imageAx'].plot([-1e5,-1e5], array([0, 1382])*p['pixelSize']*p['axScale'], 'r-', lw=1)
        self.drawCrosshair(p['lineMarker2H'], p['lineMarker2V'], p['Marker2Pos'])
        p['imageAx'].set_xlim(p['xylim'][0])
        p['imageAx'].set_ylim(p['xylim'][1])
        p['imageAx'].set_xlabel(u'\u00b5m')
        p['imageAx'].set_ylabel(u'\u00b5m')
        if not self.roiOnly_check_button.get_active():
            lc = p['partLumCorr']*p['pixelSize']*p['axScale']
            ac = p['partAtomCount']*p['pixelSize']*p['axScale']
            rectLumCorr = plt.Rectangle((lc[2], lc[0]), lc[3]-lc[2], lc[1]-lc[0], ec='red', fc='none', ls='dotted')
            rectAtomCount = plt.Rectangle((ac[2], ac[0]), ac[3]-ac[2], ac[1]-ac[0], ec='black', fc='none', ls='dashed')
            p['imageAx'].add_patch(rectLumCorr)
            p['imageAx'].add_patch(rectAtomCount)

    def initNav(self):
        home = Nav.home
        def new_home(self, *args, **kwargs):
            home(self, *args, **kwargs)
            Panoptikum_.nav_callback('home')
        Nav.home = new_home

        back = Nav.back
        def new_back(self, *args, **kwargs):
            back(self, *args, **kwargs)
            Panoptikum_.nav_callback('back')
        Nav.back = new_back

        forward = Nav.forward
        def new_forward(self, *args, **kwargs):
            forward(self, *args, **kwargs)
            Panoptikum_.nav_callback('forward')
        Nav.forward = new_forward

    def drawImageArea(self, parameters):
        p = parameters
        ImageScale = array(p['part']) * p['pixelSize'] * p['axScale']
        data = pd.Series(zeros(len(p['values2calc'])), index=p['indexLst'], dtype=float)

        # image field
        if p['imageCategory'] == 'opt. Dichte' or p['imageCategory'] == 'ROI':
            image = p['fullFrame'][p['part'][0]:p['part'][1], p['part'][2]:p['part'][3]]
            if (p['part'][1]-p['part'][0] < 800 and p['part'][3]-p['part'][2] < 800) or p['imageCategory'] == 'ROI':
                data = pd.Series(imageAnalyse.analyseShot(p['adwinfile'], p['values2calc'], image), index=p['indexLst'], dtype=float)
            imageAtomCount = p['fullFrame'][p['partAtomCount'][0]:p['partAtomCount'][1], p['partAtomCount'][2]:p['partAtomCount'][3]]
            data['n'] = imageAnalyse.analyseShot(p['adwinfile'], [p['values2calc'][0]], imageAtomCount)
            if p['ODmaxAuto'].get_active():
                p['ODmax'].set_value(imageAtomCount.max())
            if data['n'] < 1e6:
                title1 = p['imageAx'].set_title((r'N(' + p['species'] + r') = ${:,.2f}\,\times\,10^3$').format(data['n']/1e3))
                if self.MessageType == 'N_Li' and p['species'] == 'Li' or self.MessageType == 'N_Rb' and p['species'] == 'Rb':
                    self.windowMessage.set_markup(u'<span font="120">N(' + p['species'] + r')= ' + u'%.2f\u00d710\u00b3</span>' %(data['n']/1e3))
            else:
                title1 = p['imageAx'].set_title((r'N(' + p['species'] + r') = ${:,.2f}\,\times\,10^6$').format(data['n']/1e6))
                if self.MessageType == 'N_Li' and p['species'] == 'Li' or self.MessageType == 'N_Rb' and p['species'] == 'Rb':
                    self.windowMessage.set_markup(u'<span font="120">N(' + p['species'] + r')= ' + u'%.2f\u00d710\u2076</span>' %(data['n']/1e6))
        else:
            if p['imageCategory'] == 'Atomrohbild':
                image = apogee.singleImage(p['filenames'][0], p['part'])
            elif p['imageCategory'] == 'Hellrohbild':
                image = apogee.singleImage(p['filenames'][1], p['part'])
            elif p['imageCategory'] == 'Dunkelrohbild':
                image = apogee.singleImage(p['filenames'][2], p['part'])
            title1 = p['imageAx'].set_title(p['species'] + r': ' + p['imageCategory'])
            if p['ODmaxAuto'].get_active():
                p['ODmax'].set_value(image.max())
        p['imageObj'].set_data(image)
        p['imageObj'].set_extent([ImageScale[2], ImageScale[3], ImageScale[1], ImageScale[0]])
        p['imageObj'].set_clim(0, p['ODmax'].get_value())
        title1.set_fontsize(32)
        title1.set_y(1.04)
        if not data['t'] == 0 and not math.isnan(float(data['t'] * 1e6)):
            self.drawCrosshair(p['lineFitH'], p['lineFitV'], [data['posh']*p['pixelSize']*p['axScale']+ImageScale[2], data['posv']*p['pixelSize']*p['axScale']+ImageScale[0]])
        else:
            self.drawCrosshair(p['lineFitH'], p['lineFitV'], [-1e5,-1e5])

        # linescan gaussian fit
        p['linescanAx'].cla()
        p['linescanAx'].plot(linspace(ImageScale[2], ImageScale[3], size(image,1)), imageAnalyse.linescan(image))
        p['linescanAx'].plot(linspace(ImageScale[2], ImageScale[3], size(image,0)), imageAnalyse.linescan(image,dir=1), color='0.85')
        if math.isnan(float(data['t'] * 1e6)):
            title3 = p['linescanAx'].set_title(u'T(' + p['species'] + r') uneindeutig')
            p['linescanAx'].text(0.02, 0.05, u'TOF = {:.2f}ms'.format(data['abb_tof']/1e3), transform=p['linescanAx'].transAxes)
        elif data['t'] == 0:
            title3 = p['linescanAx'].set_title('')
        else:
            title3 = p['linescanAx'].set_title(u'T(' + p['species'] + u') = {:.2f}\u00b5K'.format(data['t'] * 1e6))
            p['linescanAx'].text(0.02, 0.85, u'H-Pos: {:.2f}\u00b5m'.format(data['posh']*p['pixelSize']*p['axScale']+ImageScale[2]), transform=p['linescanAx'].transAxes)
            p['linescanAx'].text(0.02, 0.7, u'H-FWHM: {:.2f}\u00b5m'.format(data['fwhmh']*p['axScale']), transform=p['linescanAx'].transAxes)
            p['linescanAx'].text(.98, 0.85, u'V-Pos: {:.2f}\u00b5m'.format(data['posv']*p['pixelSize']*p['axScale']+ImageScale[0]), transform=p['linescanAx'].transAxes, color='0.7', horizontalalignment='right')
            p['linescanAx'].text(.98, 0.7, u'V-FWHM: {:.2f}\u00b5m'.format(data['fwhmv']*p['axScale']), transform=p['linescanAx'].transAxes, color='0.7', horizontalalignment='right')
            p['linescanAx'].plot(linspace(ImageScale[2], ImageScale[3], size(image,1)), imageAnalyse.gauss(range(0,size(image.conj().transpose(),0)), data['amph'],data['posh'],imageAnalyse.fwhm(data['fwhmh'], reverse=1),data['offh']), '--')
            p['linescanAx'].plot(linspace(ImageScale[2], ImageScale[3], size(image,0)), imageAnalyse.gauss(range(0,size(image.conj().transpose(),1)), data['ampv'],data['posv'],imageAnalyse.fwhm(data['fwhmv'], reverse=1),data['offv']), '--', color='0.85')
            p['linescanAx'].plot([data['posh']*p['pixelSize']*p['axScale']+ImageScale[2], data['posh']*p['pixelSize']*p['axScale']+ImageScale[2]], p['linescanAx'].get_ylim(), 'k:')
            posVrelAxH = (data['posv']*p['pixelSize']*p['axScale'])/(ImageScale[1]-ImageScale[0])*(ImageScale[3]-ImageScale[2])+ImageScale[2]
            p['linescanAx'].plot([posVrelAxH, posVrelAxH], p['linescanAx'].get_ylim(), ':', color='0.5')
            p['linescanAx'].plot([p['Marker1Pos'][0], p['Marker1Pos'][0]], p['linescanAx'].get_ylim(), 'k-')
            posVrelAxHMarker1 = (p['Marker1Pos'][1]-ImageScale[0])/(ImageScale[1]-ImageScale[0])*(ImageScale[3]-ImageScale[2])+ImageScale[2]
            p['linescanAx'].plot([posVrelAxHMarker1, posVrelAxHMarker1], p['linescanAx'].get_ylim(), '-', color='0.85')
            p['linescanAx'].text(0.02, 0.05, u'TOF = {:.2f}ms'.format(data['abb_tof']/1e3), transform=p['linescanAx'].transAxes)
            if not data['automatisierung'] == 0:
                p['linescanAx'].text(0.98, 0.05, u'auto: {:.0f}/{:.0f} #{:.0f}'.format(data['auto_step_pos'], abs(data['auto_break']), data['auto_n_pos']), transform=p['linescanAx'].transAxes, horizontalalignment='right')
        p['linescanAx'].set_xlim((ImageScale[2], ImageScale[3]))
        p['linescanAx'].grid(True, color='0.5')
        p['linescanAx'].set_xlabel(u'\u00b5m')
        p['linescanAx'].ticklabel_format(style='sci', axis='y', scilimits=(0,3))
        title3.set_fontsize(18)
        title3.set_y(1.06)
        p['canvasObj'].draw_idle()
        return data

    def updateADwinID(self, parameters):
        p = parameters
        if not self.live_mode_check_button.get_active():
            ID = float(p['adwinID'].get_value())
        else:
            with open(p['adwinfile'] + '.id') as adwinIDfile:
                ID = float(adwinIDfile.read())
                p['adwinID'].set_range(0,1e6)
        p['adwinID'].set_value(ID)

    def updateFullFrameImage(self, parameters):
        p = parameters
        if p['imageCategory'] == 'ROI':
            p['fullFrame'] = apogee.singleImage(p['filenames'][0], self.maxImageSize)
            x, y = shape(p['fullFrame'])
            p['part'] = array([0, x, 0, y])
        else:
            p['fullFrame'] = apogee.odTriple(p['filenames'], self.maxImageSize, p['partLumCorr'], p['partAtomCount'])
        self.updateADwinID(parameters)

    def drawCrosshair(self, line1, line2, position):
        line1.set_ydata(position[1])
        line2.set_xdata(position[0])

    def on_togglebuttonLiHide_toggled(self,button):
        if button.get_active():
            self.hbox_Li.hide()
        else:
            self.hbox_Li.show()

    def on_togglebuttonRbHide_toggled(self,button):
        if button.get_active():
            self.hbox_Rb.hide()
        else:
            self.hbox_Rb.show()

    def radiobuttonRb_toggled(self, button):
        if button.get_active():
            if button.get_label() == 'opt. Dichte':
                self.RbParameters['ODmax'].set_range(0.05, 3)
                self.RbParameters['ODmax'].set_digits(2)
                self.RbParameters['ODmaxLabel'].set_text('max opt. Dichte:')
            else:
                self.RbParameters['ODmax'].set_range(0.05, 65536)
                self.RbParameters['ODmax'].set_digits(0)
                self.RbParameters['ODmaxLabel'].set_text('max Belichtung:')
            self.RbParameters['ODmaxAuto'].set_active(True)
            self.RbParameters['imageCategory'] = button.get_label()
            self.drawImageArea(self.RbParameters)

    def radiobuttonLi_toggled(self, button):
        if button.get_active():
            if button.get_label() == 'opt. Dichte':
                self.LiParameters['ODmax'].set_range(0.05, 3)
                self.LiParameters['ODmax'].set_digits(2)
                self.LiParameters['ODmaxLabel'].set_text('max opt. Dichte:')
            else:
                self.LiParameters['ODmax'].set_range(0.05, 65536)
                self.LiParameters['ODmax'].set_digits(0)
                self.LiParameters['ODmaxLabel'].set_text('max Belichtung:')
            self.LiParameters['ODmaxAuto'].set_active(True)
            self.LiParameters['imageCategory'] = button.get_label()
            self.drawImageArea(self.LiParameters)

    def radiobuttonMessage_toggled(self, button):
        if button.get_active():
            if button.get_label() == 'N(Li)':
                self.MessageType = 'N_Li'
                self.drawImageArea(self.LiParameters)
            elif button.get_label() == 'N(Rb)':
                self.MessageType = 'N_Rb'
                self.drawImageArea(self.RbParameters)

    def hscaleRbOD_adj_value_changed(self, slider):
        self.RbParameters['imageObj'].set_clim(0, slider.get_value())
        self.RbParameters['canvasObj'].draw_idle()

    def hscaleLiOD_adj_value_changed(self, slider):
        self.LiParameters['imageObj'].set_clim(0, slider.get_value())
        self.LiParameters['canvasObj'].draw_idle()

    def filechooserbutton1Rb_selection_changed(self, event):
        print self.TEST.get_filename()
        
    def checkbuttonRbAutoOD_toggled(self, button):
        if button.get_active():
            self.RbParameters['ODmax'].set_sensitive(False)
            self.drawImageArea(self.RbParameters)
        else:
            self.RbParameters['ODmax'].set_sensitive(True)

    def checkbuttonLiAutoOD_toggled(self, button):
        if button.get_active():
            self.LiParameters['ODmax'].set_sensitive(False)
            self.drawImageArea(self.LiParameters)
        else:
            self.LiParameters['ODmax'].set_sensitive(True)

    def mouse_scrolled_callback_Rb(self, event):
        self.xlimChange(event, self.RbParameters)

    def mouse_scrolled_callback_Li(self, event):
        self.xlimChange(event, self.LiParameters)

    def xlimChange(self, event, parameters):
        p = parameters
        if event.inaxes == p['imageAx']:
            if event.button == 'up':
                x0 = p['xylim'][0][0]
                x1 = p['xylim'][0][1]
                removeFromX = (x1-x0)*0.03
                mouseXrel = (event.xdata-x0) / (x1-x0)
                x0 = x0 + removeFromX*mouseXrel
                x1 = x1 - removeFromX*(-mouseXrel+1)
                y0 = p['xylim'][1][0]
                y1 = p['xylim'][1][1]
                removeFromY = (y1-y0)*0.03
                mouseYrel = (event.ydata-y0) / (y1-y0)
                y0 = y0 + removeFromY*mouseYrel
                y1 = y1 - removeFromY*(-mouseYrel+1)
            else:
                x0 = p['xylim'][0][0]
                x1 = p['xylim'][0][1]
                addToX = (x1-x0)*0.03
                mouseXrel = (event.xdata-x0) / (x1-x0)
                x0 = x0 - addToX*mouseXrel
                x1 = x1 + addToX*(-mouseXrel+1)
                y0 = p['xylim'][1][0]
                y1 = p['xylim'][1][1]
                addToY = (y1-y0)*0.03
                mouseYrel = (event.ydata-y0) / (y1-y0)
                y0 = y0 - addToY*mouseYrel
                y1 = y1 + addToY*(-mouseYrel+1)
            p['xylim'] = [(x0, x1), (y0, y1)]
            p['imageAx'].set_xlim(p['xylim'][0])
            p['imageAx'].set_ylim(p['xylim'][1])
            p['canvasObj'].draw_idle()

    def mouse_press_callback_Rb(self, event):
        return

    def mouse_press_callback_Li(self, event):
        return

    def mouse_release_callback_Rb(self, event):
        self.mouse_release_callback(event, self.RbParameters)

    def mouse_release_callback_Li(self, event):
        self.mouse_release_callback(event, self.LiParameters)

    def mouse_release_callback(self, event, parameter):
        p = parameter
        if event.button == 2:
             if event.inaxes != p['imageAx']:
                 event.xdata, event.ydata = -1e5, -1e5
             if self.shiftPressed == False:
                 p['Marker1Pos'] = (event.xdata, event.ydata)
                 self.drawCrosshair(p['lineMarker1H'], p['lineMarker1V'], p['Marker1Pos'])
             else:
                 p['Marker2Pos'] = (event.xdata, event.ydata)
                 self.drawCrosshair(p['lineMarker2H'], p['lineMarker2V'], p['Marker2Pos'])
             self.drawImageArea(p)
        else:
            self.ImageAxisChanged(p)

    def ImageAxisChanged(self, parameters, fullframe=False):
        p = parameters
        if p['xylim'] != [p['imageAx'].get_xlim(), p['imageAx'].get_ylim()]:
            if fullframe:
                if not p['imageCategory'] == 'ROI':
                    p['part'] = self.maxImageSize
                p['xylim'] = [(1, p['part'][3]*p['pixelSize']*p['axScale']-1), (p['part'][1]*p['pixelSize']*p['axScale']-1, 1)]
                p['imageAx'].set_xlim(p['xylim'][0])
                p['imageAx'].set_ylim(p['xylim'][1])
            else:
                p['xylim'] = [p['imageAx'].get_xlim(), p['imageAx'].get_ylim()]
                newX = map(int, array(p['xylim'][0])/p['pixelSize']/p['axScale'])
                newY = map(int, array(p['xylim'][1])/p['pixelSize']/p['axScale'])
                p['part'] = array([newY[::-1], newX]).flatten() + [-1,1,-1,1]
                p['part'] = p['part'].clip(1,1391)
                if p['part'][3] > 1039:
                    p['part'][3] = 1039
            data = self.drawImageArea(p)
            self.updateLast10(p, data)

    def updateLast10(self, parameters, data, isNew = False):
            p = parameters
            if p['species'] == 'Rb':
                if isNew == True:
                    p['last10'].columns = p['last10'].columns + 1
                p['last10'][0] = pd.DataFrame(data)#, index = ['n', 't', 'fwhmv', 'fwhmh', 'posv', 'posh', 'ampv', 'amph', 'offv', 'offh', 'automatisierung', 'auto_n_pos', 'auto_step_pos', 'auto_start', 'auto_step', 'auto_break', 'abb_tof'])
                if p['last10'].shape[1] > 10:
                    print p['last10']
                    del p['last10'][10]
                    print p['last10']
                df = p['last10'].copy()
                df.columns = df.columns * (-1)
                print df.loc['n']
                self.ax1Last10.cla()
                df.loc['n'].plot(ax=self.ax1Last10)
                self.canvas_RbLast10.draw_idle()

    def RectangleSelector_callback_Rb(self, eclick, erelease):
        self.RectangleSelected(self.RbParameters, eclick, erelease)

    def RectangleSelector_callback_Li(self, eclick, erelease):
        self.RectangleSelected(self.LiParameters, eclick, erelease)

    def RectangleSelected(self, parameters, eclick, erelease):
        p = parameters
        if not self.roiOnly_check_button.get_active():
            coordinates = array([eclick.ydata, erelease.ydata, eclick.xdata, erelease.xdata])
            coordinateScaled = array(map(int, coordinates/p['pixelSize']/p['axScale']))
            coordinateScaled = sort(coordinateScaled.reshape(2,2), axis=1).flatten()
            p['imageAx'].cla()
            if erelease.button == 1:
                p['partAtomCount'] = coordinateScaled
            elif erelease.button == 3:
                p['partLumCorr'] = coordinateScaled
                self.updateFullFrameImage(p)
            self.initLines(p)
            self.drawCrosshair(p['lineMarker1H'], p['lineMarker1V'], p['Marker1Pos'])
            data = self.drawImageArea(p)
            self.updateLast10(p, data)
            #print p['species'], p['partAtomCount'], p['partLumCorr']

    def file_changedRb(self, monitor, file, unknown, event):
        if event == gio.FILE_MONITOR_EVENT_CHANGES_DONE_HINT:
            self.updateFullFrameImage(self.RbParameters)
            data = self.drawImageArea(self.RbParameters)
            self.updateLast10(self.RbParameters, data, isNew = True)
            if self.www_Rb.get_active():
                self.figure_Rb.savefig(self.RbParameters['savefilename'], dpi=50)

    def file_changedLi(self, monitor, file, unknown, event):
        if event == gio.FILE_MONITOR_EVENT_CHANGES_DONE_HINT:
            self.updateFullFrameImage(self.LiParameters)
            data = self.drawImageArea(self.LiParameters)
            self.updateLast10(self.LiParameters, data, isNew = True)
            if self.www_Li.get_active():
                self.figure_Li.savefig(self.LiParameters['savefilename'], dpi=50)

    def statusbar1_text_pushed(self, message):
        ##self.statusbar.push(self.statusbar.get_context_id(message), message)
        return

    def nav_callback(self, button):
        if button == 'home':
            fullframe = True
        else:
            fullframe = False
        self.ImageAxisChanged(self.RbParameters, fullframe)
        self.ImageAxisChanged(self.LiParameters, fullframe)

    def on_key_press(self, widget, event):
        if event.state & gtk.gdk.SHIFT_MASK:
            self.shiftPressed = True

    def on_key_released(self, widget, event):
        if event.keyval == 65505:
            self.shiftPressed = False

    def open_toggle(self, event):
        self.filechooserProject.run()
        self.live_mode_check_button.set_active(False)

    def live_mode(self, widget):
        if widget.active:
            self.RbParameters['filenames'] = ['rb_recenta.fit', 'rb_recentb.fit', 'rb_recentc.fit']
            self.LiParameters['filenames'] = ['li_recenta.fit', 'li_recentb.fit', 'li_recentc.fit']
            self.monitorRbID = self.monitorRb.connect('changed', self.file_changedRb)
            self.monitorLiID = self.monitorLi.connect('changed', self.file_changedLi)
            self.sameID_check_button.set_sensitive(False)
            self.roiOnly_check_button.set_sensitive(False)
            self.roiOnly_check_button.set_active(False)
            for p in [self.RbParameters,self.LiParameters]:
                p['folder'] = 'F:\\INBOX\\Apogee\\'
                os.chdir(p['folder'])
                p['adwinfile'] = '..\\ADwin\\adwin_code_recent.acm'
                p['adwinID'].set_sensitive(False)
                self.updateFullFrameImage(p)
                p['imageAx'].cla()
                self.initLines(p)
                self.drawCrosshair(p['lineMarker1H'], p['lineMarker1V'], p['Marker1Pos'])
                data = self.drawImageArea(p)
        else:
            self.monitorRb.disconnect(self.monitorRbID)
            self.monitorLi.disconnect(self.monitorLiID)
            if self.filechooserProject.get_filename() is None:
                self.filechooserProject.run()
            self.sameID_check_button.set_sensitive(True)
            self.roiOnly_check_button.set_sensitive(True)
            self.parse_config_file()
            
    def roiOnly_toggled(self, widget):
        if widget.active:    
            self.RbParameters['filenames'] = ['rb_roi.fit']
            self.LiParameters['filenames'] = ['li_roi.fit']
            self.RbParameters['imageCategory'] = 'ROI'
            self.LiParameters['imageCategory'] = 'ROI'
            self.RbParameters['partAtomCount'] = self.maxImageSize
            self.LiParameters['partAtomCount'] = self.maxImageSize
            self.radio_buttons_imageCategory.set_sensitive(False)
        else:
            self.RbParameters['imageCategory'] = 'opt. Dichte'
            self.LiParameters['imageCategory'] = 'opt. Dichte'
            self.radio_buttons_imageCategory.set_sensitive(True)
        self.RbParameters['part'] = self.maxImageSize
        self.LiParameters['part'] = self.maxImageSize
        self.adwinID_change(self.RbParameters)
        self.adwinID_change(self.LiParameters)
    

    def project_file_clicked(self, event):
        self.filechooserProject.hide()
        
    def parse_config_file(self):
        import ConfigParser
        import codecs
        configFile = self.filechooserProject.get_filename()
        projectFolder = os.path.dirname(configFile)
        Config = ConfigParser.ConfigParser()
        try:
            Config.readfp(codecs.open(configFile, 'r', 'utf8'))
            if 'li' in Config.sections() and not self.roiOnly_check_button.get_active():
                self.LiParameters['filenames'] = Config.get('li', 'filenames').split()
                self.LiParameters['adwinfile'] = 'adwin_code.acm'
                self.LiParameters['partAtomCount'] = array(map(int, Config.get('li', 'part').split()))
                self.LiParameters['partLumCorr'] = array(map(int, Config.get('li', 'partLumCorr').split()))
                dirListLi = map(float,[name for name in os.listdir(projectFolder) if os.path.isdir(os.path.join(projectFolder,name))])
                #self.LiParameters['adwinID'].set_value(min(dirListLi))
                self.LiParameters['adwinID'].set_range(min(dirListLi),max(dirListLi))
                self.LiParameters['adwinID'].set_sensitive(True)
            else:
                print 'No Li in project folder config-file!'
            if 'rb' in Config.sections():
                self.RbParameters['filenames'] = Config.get('rb', 'filenames').split()
                self.RbParameters['adwinfile'] = 'adwin_code.acm'
                self.RbParameters['partAtomCount'] = array(map(int, Config.get('rb', 'part').split()))
                self.RbParameters['partLumCorr'] = array(map(int, Config.get('rb', 'partLumCorr').split()))
                dirListRb = map(float,[name for name in os.listdir(projectFolder) if os.path.isdir(os.path.join(projectFolder,name))])
                #self.RbParameters['adwinID'].set_value(min(dirListRb))
                self.RbParameters['adwinID'].set_range(min(dirListRb),max(dirListRb))
                self.RbParameters['adwinID'].set_sensitive(True)
            else:
                print 'No Rb in project folder config-file!'
        except:
            self.live_mode_check_button.set_active(True)
            print 'Opening of project failed!'
            
    def adwinIDRb_changed_value(self, event):
        self.adwinID_change(self.RbParameters)
        if self.sameID_check_button.get_active():
            if not self.RbParameters['adwinID'].get_value() == self.LiParameters['adwinID'].get_value():
                self.LiParameters['adwinID'].set_value(self.RbParameters['adwinID'].get_value())
        
    def adwinIDLi_changed_value(self, event):
        self.adwinID_change(self.LiParameters)
        if self.sameID_check_button.get_active():
            if not self.RbParameters['adwinID'].get_value() == self.LiParameters['adwinID'].get_value():
                self.RbParameters['adwinID'].set_value(self.LiParameters['adwinID'].get_value())
                
    def adwinID_change(self, parameters):
        p = parameters
        if not self.live_mode_check_button.get_active():
            if not self.roiOnly_check_button.get_active():
                self.parse_config_file()
            dir = os.path.dirname(self.filechooserProject.get_filename())
            p['folder'] = os.path.join(dir, str(int(p['adwinID'].get_value())))
            os.chdir(p['folder'])
            self.updateFullFrameImage(p)
            p['imageAx'].cla()
            self.initLines(p)
            self.drawCrosshair(p['lineMarker1H'], p['lineMarker1V'], p['Marker1Pos'])
            data = self.drawImageArea(p)
            
    def quit_toggle(self, event):
        self.__del__(event)

    def on_clicked_about(self, widget):
        about = gtk.AboutDialog()
        about.set_program_name('Panoptikum')
        about.set_version('0.3')
        about.set_copyright(u'\u00a9 Reinhardt A.W. Maier')
        about.set_comments('Panoptikum is a tool for analysing cold atoms images')
        about.set_website('http://www.zaehlwerk.net')
        #about.set_logo(gtk.gdk.pixbuf_new_from_file('battery.png'))
        about.run()
        about.destroy()

    def large_atomnumber(self, widget):
        if widget.active:
            self.windowMessage.show()
        else:
            self.windowMessage.hide()

    def large_atomnumber_button(self, button, *opt):
        self.large_atom_check_button.set_active(False)
        self.windowMessage.hide()
        return True

    def stay_on_top(self, widget):
        if widget.active:
            self.window.set_keep_above(True)
        else:
            self.window.set_keep_above(False)

    def __del__(self, event):
        if self.live_mode_check_button.get_active():
            with open('restoreLastSession.pickle', 'wb') as outfile:
                parameterToStore = array([self.RbParameters['part'], self.RbParameters['partAtomCount'],self.RbParameters['partLumCorr'],self.RbParameters['xylim'],self.RbParameters['Marker1Pos'],self.RbParameters['Marker2Pos'],self.LiParameters['part'], self.LiParameters['partAtomCount'],self.LiParameters['partLumCorr'],self.LiParameters['xylim'],self.LiParameters['Marker1Pos'],self.LiParameters['Marker2Pos']])
                pickle.dump(parameterToStore, outfile)
        gtk.main_quit()

Panoptikum_ = Panoptikum()
plt.close('all')
gtk.main()
