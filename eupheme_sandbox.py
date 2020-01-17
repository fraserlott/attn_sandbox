#EUPHEME Attribution Service Website back end
#Fraser Lott, Met Office Hadley Centre
#Crown Copyright.
def read_s3_data(modelname,region,startmonth,endmonth,startyear,endyear,eventyear,
                 comparator,obsset,varo):
    'Read in obs and models and process into homogeneous, useable data.'
    from fl_detn_tools import set_variables
    from region_attn import load_obs,load_process_data,load_process_eucleia_short
    from iris.cube import CubeList
    from time import time
    calcstart=time()

    obsfname,varo,varn,freq=set_variables(obsset,varo)
    #Imported settings based on observations set requested

    if endyear>eventyear:
        obs=load_obs(startyear,endyear,region,obsset,varo)
        #Observations (already anomalies)
    else:#need to include obs up to the event
        obs=load_obs(startyear,eventyear,region,obsset,varo)

    multi_all_seas_bc=CubeList()#Multi-model ensembles may be needed
    multi_nat_seas_bc=CubeList()
    modelroot='/project/detn/'
    #modelroot='/s3/informatics-eupheme/'
    for i,modeln in enumerate(modelname):
        print('Loading '+modeln)
        if endyear>=eventyear:
            modeldata=load_process_data(startyear,endyear,eventyear,
                                        startmonth,endmonth,region,obs,
                                        modelroot,modeln,varn,varo,freq,i)
        else:
            modeldata=load_process_eucleia_short(startyear,endyear,eventyear,
                                                 startmonth,endmonth,region,obs,
                                                 modelroot,modeln,varn,varo,freq,i)
        all_seas_bc,nat_seas_bc,event,obs_event,bias,all_seas,nat_seas,obs_seas=modeldata
        multi_all_seas_bc.append(all_seas_bc)
        multi_nat_seas_bc.append(nat_seas_bc)

    duration=time()-calcstart
    print('Read time %3d minutes and %2d seconds'%(duration//60,duration%60) )
    #Year pooling section will follow here eventually (paste over from region_attn).
    return multi_all_seas_bc,multi_nat_seas_bc,obs_event,obs_seas,varn


def calc_probs(multi_seas_bc,eventyear,obs_event,varn,comparator):
    from fl_detn_tools import area_mean, cubise_ensemble
    from region_attn import get_ens_prob_wrt_thresh

    distn=area_mean(cubise_ensemble(multi_seas_bc)).data[:,:-1]
    distn=distn.reshape(distn.size)#Flatten the distribution for R
    #Don't include final year
    prob=get_ens_prob_wrt_thresh(distn,eventyear,obs_event,varn,comparator)
    return distn,prob


def print_stats(all_prob,nat_prob):
    from numpy import inf
    far=1-(nat_prob/all_prob)
    if nat_prob==0:
        rr=inf
    else:
        rr=all_prob/nat_prob
    print('P_all = %5.3f, P_nat = %5.3f'%(all_prob,nat_prob))
    print('Fraction of Attributable Risk (FAR) = %5.3f'%far)
    print('Risk ratio = %5.1f'%rr)


class FLEvent:
    def __init__(self):
        "The data below will be replaced by input from the frontend, when it's finished."
        self.startyear = 1960
        self.endyear = 2012
        self.pool = False
        self.modelname = ['HadGEM3-A-N216']  # We may use multi-model ensembles later
        #self.modelname = ['CanESM2']
        "These data will be defined by widgets."
        self.region = None
        self.regionname = None
        self.startmonth = None
        self.endmonth = None
        self.eventyear = None
        self.comparator = None
        self.obsset = None
        self.varo = None
        self.event_description = None
        self.thresh_description = None
        "These data will be computed from model and obs data."
        self.multi_all_seas_bc = None
        self.multi_nat_seas_bc = None
        self.obs_event = None
        self.obs_seas = None
        self.varn = None
        self.all_distn = None
        self.all_prob = None
        self.nat_distn = None
        self.nat_prob = None

    def read_widgets(self, regionwidg, seasonwidg, yearwidg, eventwidg):
        from fl_detn_tools import giorgi_coord
        #print('Reading region')
        self.region = giorgi_coord(region_dict[regionwidg.value])
        self.regionname = regionwidg.value
        #self.region = giorgi_coord(regionwidg.value)
        #self.regionname = giorgi_name[regionwidg.value]
        #print('Reading season')
        self.startmonth, self.endmonth = season_dict[seasonwidg.value]
        #print('Reading year')
        self.eventyear = yearwidg.value
        #print('Reading event type')
        self.event_description = event_dropdown.value
        dropdown_options = {'high temperature': ('gt', 'CRUTEM4', 'temperature_anomaly'),
                            'low temperature': ('lt', 'CRUTEM4', 'temperature_anomaly'),
                            'high rainfall': ('gt', 'GPCC 2.5 degree', 'precip'),
                            'low rainfall': ('lt', 'GPCC 2.5 degree', 'precip')}
        self.comparator, self.obsset, self.varo = dropdown_options[self.event_description]
        thresh_options = {'high temperature': 'as hot or hotter',
                          'low temperature': 'as cold or colder',
                          'high rainfall': 'as wet or wetter',
                          'low rainfall': 'as dry or drier'}
        self.thresh_description = thresh_options[self.event_description]
        #print('widget read complete')

    def read_obs_and_model_data(self):
        computed = read_s3_data(self.modelname, self.region, self.startmonth, self.endmonth,
                                self.startyear, self.endyear, self.eventyear,
                                self.comparator, self.obsset, self.varo)
        self.multi_all_seas_bc, self.multi_nat_seas_bc = computed[0:2]
        self.obs_event, self.obs_seas, self.varn = computed[2:]
        self.all_distn, self.all_prob = calc_probs(self.multi_all_seas_bc, self.eventyear,
                                                   self.obs_event, self.varn, self.comparator)
        self.nat_distn, self.nat_prob = calc_probs(self.multi_nat_seas_bc, self.eventyear,
                                                   self.obs_event, self.varn, self.comparator)

    def rr(self):
        "Calculate and output risk ratio."
        if self.nat_prob == 0:
            from numpy import inf
            return inf
        else:
            return self.all_prob / self.nat_prob

# EUPHEME Attribution Service Website front end
# Fraser Lott, Met Office Hadley Centre
# Crown Copyright.

#% matplotlib inline
# Makes sure we plot to the notebook - alternatively %matplotlib notebook
#from pandas.plotting import register_matplotlib_converters

#register_matplotlib_converters()

giorgi_name = ['Antarctica',
               'Southern Australia', 'Northern Australia',
               'Amazon Basin ', 'Southern South America', 'Central America',
               'Western North America', 'Central North America', 'Eastern North America',
               'Alaska', 'Greenland',
               'Mediterranean Basin', 'Northern Europe',
               'Western Africa', 'Eastern Africa', 'Southern Africa', 'Sahara',
               'Southeast Asia', 'East Asia', 'South Asia',
               'Central Asia', 'Tibet', 'North Asia']
#from ipywidgets import Dropdown, BoundedIntText, Button, Output, Layout, IntProgress
#from IPython.display import display, Image, HTML
from bokeh.models import Select,TextInput,Div,Button,Paragraph,Slider
from bokeh.io import curdoc #show as showwidget
from bokeh.layouts import column #widgetbox
#from os import getcwd
webpath='http://www-hc/~hadlf/eupheme/'
linuxpath='/home/h02/hadlf/public_html/eupheme/'
# region_dropdown = Select(options={giorgi_name[i]: i for i in range(len(giorgi_name))},
#                            value=11, title='Giorgi region')
region_dropdown = Select(options=giorgi_name, value='Mediterranean Basin', title='Giorgi region')
region_dict={giorgi_name[i]: i for i in range(len(giorgi_name))}
#region_dropdown = Dropdown(options={giorgi_name[i]: i for i in range(len(giorgi_name))},
#                           value=11, description='Giorgi region')
season_dropdown = Select(options=['Winter (Dec-Feb)', 'Spring (Mar-May)',
                                  'Summer (Jun-Aug)', 'Autumn (Sep-Nov)'],
                           value='Summer (Jun-Aug)', title='Season')
season_dict={'Winter (Dec-Feb)': [12, 2], 'Spring (Mar-May)': [3, 5],
             'Summer (Jun-Aug)': [6, 8], 'Autumn (Sep-Nov)': [9, 11]}
# season_dropdown = Dropdown(options={'Winter (Dec-Feb)': [12, 2], 'Spring (Mar-May)': [3, 5],
#                                     'Summer (Jun-Aug)': [6, 8], 'Autumn (Sep-Nov)': [9, 11]},
#                            value=[6, 8], description='Season')
#yearbox = TextInput(value=2003, title='Year')#Is there a way to do it more like the line below?
#yearbox = BoundedIntText(value=2003, min=1990, max=2011, step=1, description='Year')
fl_event = FLEvent()  # Initialise my event class (not Andy's...)
yearbox = Slider(value=2003, start=1990, end=2011, step=1, title='Year')
event_dropdown = Select(options=['high temperature', 'low temperature',
                                 'high rainfall', 'low rainfall'],
                          value='high temperature', title='Event type')
# event_dropdown = Dropdown(options=['high temperature', 'low temperature',
#                                    'high rainfall', 'low rainfall'],
#                           value='high temperature', description='Event type')

curdoc().title='EUPHEME'
#graphs = Output(layout=Layout(height='90%'))
# eupheme_logo=Image(url="http://eupheme.eu/resources/Copy%20of%20IMG_0036.jpg")
title_font = '<font face="Arial, Helvetica, sans-serif" size="8">'
body_font = '<font face="Arial, Helvetica, sans-serif" size="3">'
page_title = Div(text='<table style="width:100%"><td><img '#width=100% height=100% '
                  + 'src="'+webpath+'eupheme_scaled.png"></td>'
                  #+ 'src="https://eupheme.eu/wp-content/uploads/2018/07/eupheme_logo-e1532426719756.jpg"></td>'
                  + '<td align="Center" valign="Center">' + title_font
                  + '<p>Attribution Sandbox</p></td>'
                  + '<td align="Center" valign="Center">'
                  + '<p>developed by</p><p><img width=200% height=200% '
                  #+ 'src="file://'+getcwd()+'/MOHC_scaled.png">'
                  + 'src="'+webpath+'MOHC_scaled.png">'
                  + '</p></td></table>')
#testpic=Div(text='<img src="http://www-hc/~hadlf/eupheme/eupheme_scaled.jpg">')
curdoc().add_root(column(page_title))#,testpic))


#timeseries_button = Button(description='How accurate are the simulations?',
#                           layout=Layout(width='50%'))
timeseries_button = Button(label='How accurate are the simulations?')

def timeseries_callback():
    from region_attn import plot_allvsnat_t_series
    #graphs.clear_output(wait=True)
    timeseries_title = Div(text=title_font + '<table style="width:200%"><td align="Center"><p>'
                            + timeseries_button.label + '</p></td></table>')
    timeseries_text = Div(text=body_font + '<table style="width:200%"><td align="Left"><p>'
                           + "We can check the set of simulations by assessing whether "
                           + "climate observations are in the same statistical range. "
                           + "By scaling the simulations’ response to the sea surface "
                           + "temperature (SST), this fit may be improved, but ideally "
                           + "the scaling factor will be around 1 (i.e. not needed)."
                           + '</p></td></table>')
    #with graphs:
    #    display(timeseries_title)
    plot_allvsnat_t_series(fl_event.multi_all_seas_bc, fl_event.multi_nat_seas_bc,
                           fl_event.obs_seas, fl_event.varo,filename=linuxpath+'timeseries.png')
    #    display(timeseries_text)
    tseriespng=Div(text='<img src="'+webpath+'timeseries.png">')
    curdoc().remove_root(bootpage)
    curdoc().add_root(column(timeseries_title,tseriespng,timeseries_text))

timeseries_button.on_click(timeseries_callback)


bootstrap_button = Button(label='How certain are we of this?')
#bootstrap_button = Button(description='How certain are we of this?',
#                          layout=Layout(width='50%'))

def bootstrap_callback():
    # Do a big bootstrap to get sample uncertainties on probabilities
    from region_attn import make_bootstrap_distn, plot_bootstrap_distn
    global bootpage
    n_bootstraps = 500
    sample_size = 525
    all_prob_distn = make_bootstrap_distn(fl_event.all_distn, fl_event.eventyear,
                                          fl_event.obs_event, fl_event.varn,
                                          fl_event.comparator, n_bootstraps, sample_size)
    nat_prob_distn = make_bootstrap_distn(fl_event.nat_distn, fl_event.eventyear,
                                          fl_event.obs_event, fl_event.varn,
                                          fl_event.comparator, n_bootstraps, sample_size)
    rr_distn = all_prob_distn / nat_prob_distn
    #graphs.clear_output(wait=True)
    boot_title = Div(text=title_font + '<table style="width:200%"><td align="Center"><p>'
                      + bootstrap_button.label + '</p></td></table>')
    boot_text = Div(text=body_font + '<table style="width:200%"><td align="Left"><p>'
                     + "By running this analysis several times, "
                     + "weighting the data differently (a technique called bootstrapping), "
                     + "the result varies a little. The graph above shows by how much."
                     + '</p></td></table>')
    #with graphs:
    #    display(boot_title)
    plot_bootstrap_distn(rr_distn, fractions=True,filename=linuxpath+'bootstrap.png')
    #    display(boot_text, timeseries_button)
    bootpng=Div(text='<img src="'+webpath+'bootstrap.png">')
    curdoc().remove_root(natpage)
    bootpage=column(boot_title,bootpng,boot_text,timeseries_button)
    curdoc().add_root(bootpage)

bootstrap_button.on_click(bootstrap_callback)


nat_distn_button = Button(label='What fraction is our fault?')
#nat_distn_button = Button(description='What fraction is our fault?',
#                          layout=Layout(width='50%'))

def nat_distn_callback():
    from region_attn import plot_allvsnat_distn
    global natpage
    #graphs.clear_output(wait=True)
    nat_title = Div(text=title_font + '<table style="width:200%"><td align="Center"><p>'
                     + nat_distn_button.label + '</p></td></table>')
    nat_text = Div(text=body_font + '<table style="width:200%"><td align="Left"><p>'
                    + 'We examine climate simulations with and without human influences '
                    + 'to get two alternative distributions. The ratio of the areas beyond '
                    + "the threshold shows us how much the event’s probability has changed."
                    + '</p></td></table>')
    #with graphs:
    #    display(nat_title)
    plot_allvsnat_distn(fl_event.all_distn, fl_event.obs_event, fl_event.nat_distn,
                        fl_event.varo, fl_event.varn, fl_event.comparator,filename=linuxpath+'nat_distn.png')
    print_stats(fl_event.all_prob, fl_event.nat_prob)
    # Replace the above with easier data output
    #    display(nat_text, bootstrap_button)  # Button to display for bootstrap distn
    natpng=Div(text='<img src="'+webpath+'nat_distn.png">')
    curdoc().remove_root(allpage)
    natpage=column(nat_title,natpng,nat_text,bootstrap_button)
    curdoc().add_root(natpage)

nat_distn_button.on_click(nat_distn_callback)


all_distn_button = Button(label='How do we calculate this?')
#all_distn_button = Button(description='How do we calculate this?',
#                          layout=Layout(width='50%'))
def all_distn_callback():
    from region_attn import plot_allvsnat_distn
    global allpage
    #graphs.clear_output(wait=True)
    all_title = Div(text=title_font + '<table style="width:200%"><td align="Center"><p>'
                     + all_distn_button.label + '</p></td></table>')
    all_text = Div(text=body_font + '<table style="width:200%"><td align="Left"><p>'
                    + "This is the distribution of all summer temperatures in this region. "#TODO:Parse data for text
                    + "We use this year’s " + fl_event.event_description + " as a threshold, "
                    + "and we’re interested in any event that could have been "
                    + fl_event.thresh_description + ".</p></td></table>")
    #with graphs:
        #display(all_title)
    plot_allvsnat_distn(fl_event.all_distn, fl_event.obs_event, None,
                        fl_event.varo, fl_event.varn, fl_event.comparator,filename=linuxpath+'all_distn.png')
        #display(all_text, nat_distn_button)  # Button to display the next graph
    allpng=Div(text='<img src="'+webpath+'all_distn.png">')
    curdoc().remove_root(barpage)
    allpage=column(all_title,allpng,all_text,nat_distn_button)
    curdoc().add_root(allpage)

all_distn_button.on_click(all_distn_callback)

# # Set up text buffer to allow us to monitor the output for the progress bar
# from io import StringIO
# from contextlib import redirect_stdout
#
# buffer = StringIO()
#
# progressbar = IntProgress(value=0, min=0, max=32, step=1,
#                           description='Progress:', layout=Layout(color='#e66f0e'))
# progressbar.style.bar_color = '#e66f0e'
#
#
# # Max should be number of print lines produced by read_s3_data.
# # Without KubeCluster, should be 32 for long runs only, 1077 if we have Ext runs
# # With KubeCluster, should be 6.
# def poll_buffer():
#     'Update the progress bar according to the length of the output buffer.'
#     from time import sleep
#     bufferlen = 0
#     buffer.flush()
#     timetaken = 0
#     while bufferlen < progressbar.max and timetaken < 13 * 60:
#         progressbar.value = bufferlen  # Update progress bar
#         sleep(1)  # second
#         timetaken += 1
#         bufferlen = len(buffer.getvalue().splitlines())  # Count lines of print output
#

calc_button = Button(label='Calculate')


def read_dropdown_callback():
    from calendar import month_name
    global barpage
    fl_event.read_widgets(region_dropdown, season_dropdown, yearbox, event_dropdown)
    # graphs.clear_output()
    # with graphs:
    #     print('Analysing ' + fl_event.regionname, fl_event.region)
    #     print(fl_event.comparator + ' %4d ' % fl_event.eventyear
    #           + month_name[fl_event.startmonth] + ' to '
    #           + month_name[fl_event.endmonth] + ' ' + fl_event.obsset + ' ' + fl_event.varo)
    #     display(progressbar)
    print('Analysing ' + fl_event.regionname,fl_event.region)
    print(fl_event.comparator + ' %4d ' % fl_event.eventyear
          + month_name[fl_event.startmonth] + ' to '
          + month_name[fl_event.endmonth] + ' ' + fl_event.obsset + ' ' + fl_event.varo)
    #showwidget(widgetbox([progressbar]))
    # Now use the text output of fl_event.read_obs_and_model_data() to monitor progress
    from threading import Thread
    # Parallel thread monitoring the buffer and updating the progress bar
    #pollthread = Thread(target=poll_buffer)
    #pollthread.start()
    print('Reading data')
    #with buffer, redirect_stdout(buffer):
    fl_event.read_obs_and_model_data()

    from fl_detn_tools import round2sigfigs
    rr = fl_event.rr()

    if 0.667 < rr < 1.5:
        rrmoreless = 'about as'  # or use percentages. Also guage with bootstrap.
    elif rr > 1:
        rrmoreless = 'around %gx more' % round2sigfigs(rr, 1)
    elif rr < 1:
        rrmoreless = 'around %gx less' % round2sigfigs(1. / rr, 1)
    rr_text = Div(text=body_font + '<table style="width:200%"><td align="Center"><p>'
                   + "This event is " + rrmoreless + " likely with human influence "
                   + 'on the climate.</p></td></table>')
    #graphs.clear_output(wait=True)
    curdoc().remove_root(selection_widgets)
    # with graphs:
    #     from event_def import Event as ACEvent  # Andy's event class
    #     from qfigs import barchart_p1p0  # Andy's simple barchart routine
    #     from matplotlib.pyplot import show
    #     try:
    #         ac_event = ACEvent([fl_event.varn, 'season', fl_event.regionname, fl_event.obs_event],
    #                            fl_event.nat_prob, fl_event.all_prob)  # Include period logic later
    #         fig, ax = barchart_p1p0([ac_event], bounds=False)  # Include adjacent regions?
    #         show(block=False)
    #
    #         display(rr_text, all_distn_button)
    #     except TypeError:
    #         print("Couldn't plot, data is missing")

    from event_def import Event as ACEvent  # Andy's event class
    from qfigs import barchart_p1p0  # Andy's simple barchart routine
    from matplotlib.pyplot import close
    try:
        ac_event = ACEvent([fl_event.varn, 'season', fl_event.regionname, fl_event.obs_event],
                           fl_event.nat_prob, fl_event.all_prob)  # Include period logic later
        fig, ax = barchart_p1p0([ac_event], bounds=False)  # Include adjacent regions?
        fig.savefig(linuxpath+"sandbox_barchart.png")
        close(fig)
        #barpng=Div(text='<img width=80% height=80% src="file://'+getcwd()+'/sandbox_barchart.png">')
        barpng=Div(text='<img src="'+webpath+'sandbox_barchart.png">')
        barpage=column(barpng,rr_text, all_distn_button)
        curdoc().add_root(barpage)
    except TypeError:
            print("Couldn't plot, data is missing")


calc_button.on_click(read_dropdown_callback)
selection_widgets=column(region_dropdown, season_dropdown, yearbox, event_dropdown, calc_button)
curdoc().add_root(selection_widgets)


# display(page_title,
#        region_dropdown,season_dropdown,yearbox,event_dropdown,
#        calc_button,graphs)#This version keeps dropdowns and button above graphs
#display(page_title, graphs)  # This version has them be replaced by the graphs
# with graphs:
#     display(region_dropdown, season_dropdown, yearbox, event_dropdown, calc_button)
