'''
Created on May 29, 2014

@author: Fraser Lott, Met Office Hadley Centre. Crown Copyright.
'''
import numpy as np
import iris.analysis as ia

domain_dict={'psl':'A','ps':'A','tas':'A','pr':'A','zg':'A','ts':'L',
             'mrsos':'L','mrso':'L','mrros':'L','mrrob':'E',
             'sfcWind':'A','ua':'A','va':'A','wap':'A',
             'hfls':'A','rsus':'A','rsds':'A','hurs':'A','huss':'A'}#prefix for monthly domain



def repair_coord_system(model):
    'Repair issue with Iris 1.7 regridding.'
    from iris.coord_systems import GeogCS
    from iris.fileformats.pp import EARTH_RADIUS
    um_cs=GeogCS(EARTH_RADIUS)
    for coord_name in ('latitude', 'longitude'):
        coord = model.coord(coord_name)
        #assert coord.coord_system is None #Didn't work as of testing on 2017-08-09
        coord.coord_system = um_cs    


def giorgi_coord(region):
    'Coordinates of Giorgi regions (West, South, East, North)'
    coord_regs=[[0,-90,360,-60],#'Antarctica'         0
               [110,-50,180,-28],#'Southern Australia'    1
               [110,-28,155,-10],#'Northern Australia'    2
               [-80,-20,-35,10],#'Amazon Basin '        3
               [-75,-55,-40,-20],#'Southern South America'    4
               [-115,10,-85,30],#'Central America'    5
               [-130,30,-105,60],#'Western North America'    6
               [-105,30,-85,50],#'Central North America'     7
               [-85,25,-60,50],#'Eastern North America'    8
               [-170,60,-105,70],#'Alaska'         9
               [-105,50,-10,85],#'Greenland'        10
               [-10,30,40,50],#'Mediterranean Basin'    11
               [-10,50,40,75],#'Northern Europe'     12
               [-20,-10,20,20],#'Western Africa'        13
               [20,-10,50,20],#'Eastern Africa'        14
               [-10,-35,50,-10],#'Southern Africa'    15
               [-20,20,65,30],#'Sahara'            16
               [95,-10,155,20],#'Southeast Asia'        17
               [100,20,145,50],#'East Asia'        18
               [65,5,100,30],#'South Asia'        19
               [40,30,75,50],#'Central Asia'        20
               [75,30,100,50],#'Tibet'            21
               [40,50,180,70],# 'North Asia'        22
	       [300,-35,310,-22.5]]#SESA
    return coord_regs[region]


def trieste_regions_coord(region):
    'Coordinates of SESA and SAf regions (West, South, East, North)'
    coord_regs=[[-60,-35,-50,-22.5],#SESA
                [-10,-35,52,-12]]#SAF
    return coord_regs[region]


def get_srex_mask(region):
    "Given an SREX region number, it imports the numpy masked array that covers it."
    "(i.e. True is where we want to mask the grid box because it's not in the region.)"
    "You may find this both slow and memory-intensive, but your other option"
    "is using Ascend, and the licence isn't available outside the Met Office."
    #We could also include loading the region name.
    from iris import load_cube
    cube=load_cube('/data/cr1/aciav/srex/*%02d_mask.pp'%region)
    #repair_coord_system(cube)
    return cube


def set_variables(obsset,varo,daily=False):
    'Deduce from the obs set (and variable in the case of E-OBS) where the file is, the frequency and what the CF-compliant name is.'
    from os import getcwd
    if getcwd().startswith('/home/jovyan') or getcwd().startswith('/home/ec2-user'):#then we're running on Pangeo/AWS and the files will be on S3
        obsfname_dict={'E-OBS': '/s3/informatics-eupheme/obs/%s_0.50deg_reg_v17.0.nc'%varo,
                       'GPCC 2.5 degree': '/s3/informatics-eupheme/obs/gpcc/full_data_monthly_v2018_25.nc',
                       'CRUTEM4': '/s3/informatics-eupheme/obs/CRUTEM.4.6.0.0.anomalies.nc'}
    else:#we're on the Met Office Linux system
        obsfname_dict={'GHCNDEX_Rx1day': '/data/local/hadlf/GHCNDEX/GHCND_Rx1day_1951-2015_RegularGrid_global_2.5x2.5deg_LSmask.nc',
                      'GPCC': '/data/local/hadlf/precip.mon.total.v6.nc',
                      'GPCC 2.5 degree': '/data/cr1/hadlf/gpcc/full_data_monthly_v2018_25.nc',
                      'GPCC 1.0 degree': '/data/cr1/hadlf/gpcc/full_data_monthly_v2018_10.nc',
                      'aphrodite': '/data/cr1/hadlf/aphrodite/*.nc',#'/data/local/hadlf/obs/aphrodite/*.nc',
                      'TAMSAT': '/data/local/hadlf/data_els028/tarcat/monthly/rfe*.nc',
                      'E-OBS': '/data/cr1/hadlf/eobs/%s_0.50deg_reg_1950-2018.nc'%varo,#_v17.0.nc'%varo,
                      'SWBM': '/data/local/hadlf/obs/SWBM_soilmoisture_EOBS_1984-2013_monthly.nc',
                      'CRUTEM4': '~/detn/obs/CRUTEM.4.6.0.0.anomalies.nc'}   #Identified by the obs variable name, give file locations
    obsfname=obsfname_dict[obsset]
    if obsset=='E-OBS':
        varn_dict_eobs={'pp':'psl','rr':'pr','tg':'tas'}
        varn=varn_dict_eobs[varo]
        if daily:
            freq='day'
        else:
            freq='Amon'
            obsfname=obsfname[0:-3]+'mon.nc'
    elif obsset.startswith('GPCC'):
        varo='precip'
        varn='pr'
        freq='Amon'
    else:
        varo_dict={'CRUTEM4':'temperature_anomaly','aphrodite':'precip','TAMSAT':'rfe','GHCNDEX_Rx1day':'Ann','SWBM':'soil_moisture'}
        varo=varo_dict[obsset]
        varn_dict={'CRUTEM4':'tas','aphrodite':'pr','TAMSAT':'pr','GHCNDEX_Rx1day':'pr','SWBM':'mrso'}
        varn=varn_dict[obsset]
        freq_dict={'CRUTEM4':'Amon','aphrodite':'Amon','TAMSAT':'Amon','GHCNDEX_Rx1day':'Rx1day','SWBM':'Lmon'}
        freq=freq_dict[obsset]
    return obsfname,varo,varn,freq


def get_srex_vertices(region):
    "Searches Andy's SREX file for the numbered region and outputs a list"
    "with the title and code followed by tuples defining the vertices of the region."
    #At some point, this might be useful to update to deal with a list of regions.
    with open('/data/cr1/aciav/srex/SREX_regs.dat') as fileobj:
        for line in fileobj:
            meta_and_coords=line[0:-1].split(' & ')#Separate out labels, remove endofline
            if float(meta_and_coords[2])==region:
                title,code,_,_=meta_and_coords
                coordstr=meta_and_coords[3].split('(')#Now broken into individual coord pair strings
                coordlist=[title,code]#Coords will be appended to this
                for coord in coordstr[1:]:
                    coordnsew=coord[0:-2].split(', ')#Now down to numbers N, S, E, W
                    coordlist.append( ( nsew2signedfloat(coordnsew[1]), nsew2signedfloat(coordnsew[0]) ) )#Reversed the 1 and 0 (lat and lon)
                return coordlist#For use with Ascend, this should be a list of tuples


def get_srex_shape(region):
    'Uses Ascend to create a shape object of the polygon of an SREX region. Will not work outside the Met Office.'
    from ascend import shape
    regioncoords=get_srex_vertices(region)
    return shape.create(regioncoords[2:], {'name': regioncoords[0], 'code':regioncoords[1], 'number':str(region)}, 'Polygon')    


def get_rowell_kenya_somalia():
    from ascend import shape
    #regioncoords=[(35,3),(40,3),(44,7),(44,10.5),(50,10.5),(50,5.5),(38,-5.5),(35,-5.5)]#THESE ARE BACKWARDS!
    regioncoords=[(3,35),(3,40),(7,44),(10.5,44),(10.5,50),(5.5,50),(-5.5,38),(-5.5,35)]
    return shape.create(regioncoords,{'name':'Kenya-Somalia','number':7},'Polygon')


def get_rectangular_region_shape(region):
    'Uses Ascend to create a shape object of a custom rectangle. Will not work outside the Met Office.'
    'Not currently functioning correctly...'
    from ascend import shape
    #eg. region=[-10,35,25,45]#southern Europe = -10 to 25E, 35 to 45N
    regioncoords=[(region[1],region[0]),(region[3],region[0]),(region[3],region[2]),(region[1],region[2])]
    return shape.create(regioncoords, {'name': 'Custom rectangular region'}, 'Polygon')    


def nsew2signedfloat(nsewstr):
    'Assumes the last character of a coord indicates North, South, East or West.'
    'Convert the other characters to a float, then for S or W, negate it.'
    ordinate=float(nsewstr[0:-1])#Trim off compass point and convert to number
    if (nsewstr[-1]=='S')|(nsewstr[-1]=='W'):
        ordinate*=-1#Convert compass point to positive or negative ordinate
    return ordinate


def round2sigfigs(nparr,sigfigs=1):
    'Round a numpy array to a number of significant figures (numpy can only do decimal places).'
    #How many digits to we have? Simplest way is to take base-10 log, add 1 and round down.
    digits=np.floor(np.log10(np.abs(nparr))+1)
    return np.round(nparr,np.int32(sigfigs)-np.int32(digits))


def area_mean(cube):
    '''Produces an area weighted mean over all latitudes and longitudes.'''
    from iris.analysis.cartography import area_weights
    return cube.collapsed(['latitude','longitude'],ia.MEAN,weights=area_weights(cube))


def ens_area_mean(cube_ens):
    from iris.cube import CubeList
    area_avg_ens=CubeList()
    for member in cube_ens:
        area_avg_ens.append(area_mean(member))
    return area_avg_ens


def gaussian_distn(x,mean,stdev):
    'Surely this exists somewhere else?'
    a=x-mean
    return np.exp(-a*a/2/stdev/stdev)/stdev/np.sqrt(2*np.pi)


def gaussian_lt_thresh(distn,thresh):
    from math import erf,sqrt
    stdev=distn.std()
    mean=distn.mean()
    #Use cumulative normal distribution up to the threshold. Simply 1 minus this for prob beyond the threshold.
    #For rainfall this should be a gamma distn. Arguably this would also be better with extreme value distns.
    prob=0.5*( 1+erf( (thresh-mean)
                      /sqrt(2*stdev*stdev) ))
    return prob#1-prob to give prob>obs (for heatwaves)


def r_fit_gpd(distn_data,thresh,npy=1):#The threshold is for the variable (e.g. T) not the probability
    "Fit data past a given threshold to a generalised Pareto distribution."
    "Use the R routine, because the SciPy fitting routine doesn't seem to work."
    #We should temporarily kill stdout, because the R routine is far too verbose.
#    import sys           
#    stdout_old=sys.stdout
#    sys.stdout=open('gpd.log','w')#Keep the print output in case it might be useful
    
    from rpy2.robjects.packages import importr
    from rpy2.robjects import FloatVector
    ismev=importr('ismev')
    
    distn_mean=distn_data.mean()
    if distn_mean<=thresh:#This is a high threshold
        distn4r=FloatVector(distn_data)
    else:
        distn4r=FloatVector(distn_mean-distn_data)
    params=ismev.gpd_fit(distn4r,thresh,npy,show=False)#npy is number of obs per block. Can we bring this in too?
    scale,shape=params[9]#This should give $mle, the maximum likelihood estimate of the parameters ([scale,shape]?)
    
#    sys.stdout=stdout_old#Now turn printing back on!
    return scale,shape


def r_fit_gev(distn_data):
    from rpy2.robjects.packages import importr
    from rpy2.robjects import FloatVector
    ismev=importr('ismev')
    distn4r=FloatVector(distn_data)#Assume for now it's positive
    
    from StdOutWrapper import StdOutWrapper#Stephen Haddad's wrapper for suppressing print output. 
    with open('gev.log','w') as my_file:
        with StdOutWrapper(my_file):#Here as demo, can use show=False with this function
            params=ismev.gev_fit(distn4r)
    location,shape,scale=params[6]
    return location,shape,scale


def r_fit_covariate_gev(distn_data,thresh,covariate=None,event_value=None,obs_prob=None,event_covariate=None,type='PP'):
    from rpy2.robjects.packages import importr
    from rpy2.robjects import FloatVector,DataFrame,BoolVector
    #rstats = importr('stats', robject_translations = {"format_perc": "format_perc_R"})
    extRemes=importr('extRemes')
    distn4r=FloatVector(distn_data)
    #distn_rdict=DataFrame({'z':distn4r})
    thresh4r=FloatVector([thresh])
    #if str(type(covariate))=="<type 'bool'>":
    if covariate==None:#Natural forcings, no covariate
        fit=extRemes.fevd(distn4r,threshold=thresh4r,type=type,time_units='years')
        "type defaults to Poisson Process as per Jeon et al 2016. Use type='DF' for GEV, type='GP' for generalised Pareto distn."
    else:#All forcings (or obs) uses global mean as a covariate on the location parameter
        global4r=FloatVector(covariate)
        global_rdict=DataFrame({'x':global4r})
        fit=extRemes.fevd(distn4r,global_rdict,threshold=thresh4r,location_fun='~x',type=type,time_units='years')#,time_units='m/year')
        #Check this works right for low extremes
    
    if obs_prob==None:#then calculate the probability wrt the observed value. (GEV concerns block maxima, eg. tasmax or tasmin, not tas)
        lowbool=BoolVector([thresh<0.5])
        event4r=FloatVector([event_value])
        if event_covariate==None:
            prob_wrt_obs=extRemes.pextRemes(fit,event4r,lower_tail=lowbool)
        else:
            eventcov4r=FloatVector([event_covariate])
            prob_wrt_obs=extRemes.pextRemes(fit,event4r,qcov=eventcov4r,lower_tail=lowbool)
        return np.array(prob_wrt_obs),fit
    elif obs_prob>0:#get return level (e.g. T or pr) for a given return period (we must convert this from prob as that's what we have).
        returntime=FloatVector([1./obs_prob])#assuming that was the prob of happening on a given year
        if event_covariate==None:
            level=extRemes.return_level(fit,returntime)
        else:
            eventcov4r=FloatVector([event_covariate])
            level=extRemes.return_level(fit,returntime,qcov=eventcov4r)
        ##Better to pass fit out and apply to return_level() or period later? We might need both. Or have the option of which is passed in?
        return np.array(level),fit
    else:
        return None,fit


def r_fit_gamma(distn_data,distn_name='gamma'):
    '''Given a distribution of data, uses R to fit, by default, a Gamma distribution.
    This works a LOT more reliably than the home-made version or the Scipy version.
    This is presumably because of a combination of a lack of documentation and understanding.'''
    from rpy2.robjects.packages import importr
    from rpy2.robjects import FloatVector
    rmass=importr('MASS')
    params=rmass.fitdistr(FloatVector(distn_data),distn_name)
    shape,rate=params[0]
    return shape,rate


def r_gamma_distn(xaxis,shape,rate):
    '''Given a shape and rate parameter and the points on the x axis, computes the Gamma PDF.
    This works a LOT more reliably than the home-made version or the Scipy version.
    This is presumably because of a combination of a lack of documentation and understanding.'''
    from rpy2.robjects.packages import importr
    from rpy2.robjects import FloatVector
    rstats=importr('stats')
    probden=rstats.dgamma(FloatVector(xaxis),shape,rate=rate)
    return np.array(probden)


def r_gamma_lt_thresh(shape,rate,thresh):
    '''Get the probability of exceeding the threshold in a gamma distribution
    of a given shape and rate. More certain this will work than the Scipy version.'''
    from rpy2.robjects.packages import importr
    rstats=importr('stats')
    prob=rstats.pgamma(thresh,shape=shape,rate=rate)
    return prob


def reciprocal_label(tick):
    '''From a matplotlib.axis.XTick object or a float,
    return a string representation x, or as 1/x for numbers less than 1.'''
    from matplotlib.axis import XTick
    if type(tick)==XTick:
        tickval=tick.get_loc()
    else:
        tickval=tick
        
    if tickval==0:
        return '-inf'
    elif tickval<1:
        return '1/%g'%(1./tickval)
    else:
        return '%g'%tickval


def label_dates(cube):
    'Adds auxiliary time coordinates "day_of_year", "month_number" and "year"'
    from iris.coord_categorisation import add_year,add_month_number,add_day_of_year
    add_year(cube, 'time', name='year')
    add_month_number(cube, 'time', name='month_number')
    add_day_of_year(cube, 'time', name='day_of_year')


def remove_datelabels(cube):
    'Remove auxiliary time coordinates "day_of_year", "month_number" and "year"'
    cube.remove_coord("day_of_year")
    cube.remove_coord("month_number")
    cube.remove_coord("year")
    

def datelabel_callback(cube,field,fname):
    if cube.coords('year')==[]:
        label_dates(cube)


def filter_months_old(cube_in,startmonth,endmonth):
    '''Returns the input cube with only the desired months of the year present.
    Could also be done with iris.Constraint().'''
    month=cube_in.coord('month_number').points
    if startmonth<=endmonth:#Doesn't cross new year
        selected = np.logical_and(startmonth<=month,month<=endmonth)
        cube_out=cube_in[selected]#This assumes that the time coordinate is the first axis of the cube.
    else:#Does cross new year
        nextyear = month<=endmonth
        selected = np.logical_or(startmonth<=month,nextyear)
        cube_copy=cube_in.copy()#So we don't overwrite the months on the original data
        cube_copy.coord('year').points[nextyear]-=1#Shifts the new year to after the end of the season
        cube_out=cube_copy[selected]#This assumes that the time coordinate is the first axis of the cube.
    return cube_out


def filter_months(cube_in,startmonth,endmonth,yearshift=True):
    '''Returns the input cube with only the desired months of the year present.
    This uses iris.Constraint(). Caveat emptor: If your season straddles the new year,
    this code (and the old version) will rewrite the 'year' coord to be the same year
    throughout the season. This is the same as year of the start month.
    To avoid this functionality, set yearshift=False.'''
    from iris import Constraint
    if startmonth<=endmonth:#Doesn't cross new year
        inseason=Constraint(coord_values={'month_number': lambda month: startmonth<=month<=endmonth})
        return cube_in.extract(inseason)
    else:
        #inseason=Constraint(coord_values={'month_number': lambda month: np.logical_and(endmonth<=month,month<=startmonth)})
        inseason=Constraint(coord_values={'month_number': lambda month: np.logical_or(startmonth<=month,month<=endmonth)})
        cube_copy=cube_in.copy()#So we don't overwrite the months on the original data
        if yearshift:
            nextyear = cube_in.coord('month_number').points<=endmonth
            cube_copy.coord('year').points[nextyear]-=1#Shifts the new year to after the end of the season
        return cube_copy.extract(inseason)


def prep_for_merge(cube):#,member_axis='realization'):
    from iris.coords import DimCoord, AuxCoord
    from iris.exceptions import CoordinateNotFoundError
    for member_axis in ['realization','physics_version']:
        if member_axis in cube.attributes:
            ens=cube.attributes.pop(member_axis)
            cube.add_aux_coord(DimCoord(np.int32(ens),long_name=member_axis))#'realization'))
    for metatitle in ['history','creation_date','branch_time','tracking_id','mo_runid',
                      'CCCma_runid','CCCma_parent_runid','parent_experiment','parent_experiment_rip',
                      'cmor_version','table_id','forcing','contact','source',
                      'branch_time_YMDH','associated_files','gfdl_experiment_name', 'title',
                      'cesm_casename','cesm_repotag', 'cesm_compset','processed_by','processing_code_information']:#member_axis,
        #if cube.attributes.has_key(metatitle):#Deprecated in Python3
        if metatitle in cube.attributes:
            metadata=cube.attributes.pop(metatitle)
            #cube.add_aux_coord(AuxCoord(metadata,long_name=metatitle))#Temporarily remove this if causing problems with merging
    for metatitle in ['comment','references']:#Just delete these ones
        #if cube.attributes.has_key(metatitle):#Deprecated in Python3
        if metatitle in cube.attributes:
            metadata=cube.attributes.pop(metatitle)
    try:
        cube.remove_coord('height')
    except CoordinateNotFoundError:
        #print('height not present')
        pass


def prep_for_merge_callback(cube,field,fname):
    'Callback version of the above.'
    prep_for_merge(cube)



def cubise_ensemble(cubelist,member_axis='realization',leavelist=False):
    'Turns a cubelist ensemble into a cube, with "realization" as the ensemble axis'
    n_memb=np.empty(len(cubelist))
    for i,cube in enumerate(cubelist):
        prep_for_merge(cubelist[i])#,member_axis)
        n_memb[i]=cube.coord(member_axis).points.size#Note the number of members for each cube
        #This should be 1 in most cases, unless we're merging multiple ensembles (e.g. multimodel mode)
    if leavelist:
        return cubelist.merge()
    elif (n_memb>1).sum() > 1:#More than one cube has more than one ensemble member
        for cube in cubelist:
            cube.remove_coord('day_of_year')
        from iris.cube import CubeList
        cubelist_out=CubeList([cubelist[0]])
        #print(cubelist[0].coord(member_axis).points)
        for i in range(1,len(cubelist)):
            #Make member numbers continue from the end of the ensemble to which we are attaching it.
            cube_out=cubelist[i].copy()
            #print('Appending to member axis')
            cube_out.coord(member_axis).points = cubelist[i].coord(member_axis).points + cubelist_out[i-1].coord(member_axis).points.max()
            #print(cube_out.coord(member_axis).points)
            cubelist_out.append(cube_out)
        #global concat_return
        #concat_return = cubelist_out
        return cubelist_out.concatenate_cube()
    else:
        return cubelist.merge_cube()


def calc_climatology(cube_in,startyear,endyear):
    '''Works out climatology for a monthly series cube, between startyear and endyear.
    This version uses constraints. Should work for other timescales too.'''
    from iris import Constraint
    #inclim=Constraint(coord_values={'year': lambda year: np.logical_and(startyear<=year,year<=endyear)})
    inclim = Constraint(year = lambda y: (y >= startyear) & (y <= endyear))#Seems to work when the above doesn't
    cube_inclim=cube_in.extract(inclim)
    #from pdb import set_trace
    #set_trace()
    if len(cube_inclim.coord('time').points) > (endyear-startyear+1)*12:#then these are dailies, not monthlies
        return cube_inclim.aggregated_by('day_of_year', ia.MEAN)#so do daily climatology
    else:#do monthly climatology
        return cube_inclim.aggregated_by('month_number', ia.MEAN)


def constrain_date(startmonth,endmonth,startyear,endyear):
    'Iris 2.1 version to get cube constraint from the start and end month and year. No auxiliary axes required.'
    from iris import Constraint
    if startyear==endyear:
        #return Constraint(month_number=lambda month: startmonth<=month<=endmonth)#Requires labelled dates
        return Constraint(time=lambda cell: startmonth<=cell.point.month<=endmonth)
    elif endyear>startyear:
        return Constraint(time=lambda cell: np.logical_or(np.logical_and(cell.point.year>startyear,cell.point.year<endyear),
                                                          np.logical_or(np.logical_and(cell.point.year==startyear, cell.point.month>=startmonth),
                                                                        np.logical_and(cell.point.year==endyear, cell.point.month<=endmonth)      )) )
        #Perhaps easier to understand below:                                                   
        #(year>startyear and year<endyear) or (year==startyear and month_number>=startmonth) or (year==endyear and monthnumber<=endmonth)
    else:
        raise ValueError('Your start year is after your end year...')



def anomalise(cube_in,startyear,endyear,clim=None):
    '''Returns the anomaly of a cube (based on T method) for a climatology between startyear and endyear.
    This is unlikely to work well if there is missing time series data.'''
    year=cube_in.coord('year').points
    cubemaxyear=year.max()
    cubeminyear=year.min()
    n_years=cubemaxyear-cubeminyear+1
    cube_out=cube_in.copy()
    from iris.cube import Cube
    #from dask import delayed
    if type(clim)==Cube or str(type(clim))=="<class 'dask.delayed.DelayedLeaf'>":
        climatology=clim
    else:#if clim==None:
        climatology=calc_climatology(cube_in,startyear,endyear)

    if len(climatology.coord('time').points)>12:#then it's a daily climatology
        cubecalendar=cube_in.coord('time').units.calendar
        if cubecalendar=='gregorian' or cubecalendar=='standard':#then deal with leap years
            from calendar import isleap
            npleap=np.vectorize(isleap)
            leapyear=npleap(year)#Just the points in years that are leap years
            n_leaps=npleap(np.arange(cubeminyear,cubemaxyear+1)).sum()
            #Now here's a big assumption - that the time data is in the right order...
            cube_out.data[leapyear]-=np.tile(climatology.data, (n_leaps,1,1) )
            #cube_out[np.logical_not(leapyear)].data-=np.tile(np.delete(climatology.data,59,0), (n_years-n_leaps,1,1) )
            #Removed 29th Feb from the climatology. Or should this be 31st Dec? Try that below
            cube_out.data[np.logical_not(leapyear)]-=np.tile( climatology.data[0:365], (n_years-n_leaps,1,1) )
        else:#don't bother with leap years - this probably needs refining!!
            cube_out.data-=np.tile(climatology.data, (n_years,1,1) )#Assumes time is always the first dimension of a cube!
    else:#monthly climatology
        #cube_out.data-=np.tile(climatology.data, (n_years,1,1) )#Assumes time is always the first dimension of a cube!
        #Check which dimension is time
        repeat=np.ones(len(cube_in.shape),dtype='int32')#Default to 1 (i.e. no repeat) in each of the cube's dimensions
        timedim=cube_in.coord_dims('time')[0]
        repeat[timedim]=n_years#Repeat for the number of years in the time dimension
        repeated_clim=np.tile(climatology.data, repeat)

        #Check if it starts in January
        janstart=np.where(cube_in.coord('month_number').points==1)[0][0]#The first January element
        selection=np.arange((12-janstart)%12, repeated_clim.shape[timedim])#Trim the climatology down for starting in e.g. December
        #print(cube_in.shape)
        #print(cube_in.coord('realization').points)
        #print(repeated_clim.take(selection,axis=timedim).shape)
        #If we need to go into the next year too, we'll need another repeat on the tile.
        cube_out.data-=repeated_clim.take(selection,axis=timedim)#could also do this with nparray.compress if selection was boolean
    return cube_out


# def harmonise_time(cubelist):
#     "Convert all cubes' time units ('time since...') to that of the first cube in the list."
#     basetimeunit=cubelist[0].coord('time').units
#     if cubelist[1].coord('time').units != basetimeunit:#Assume the rest aren't right either
#         cubelist[0].coord('time').points=cubelist[0].coord('time').points.astype('float64')
#         #Because the rest will now be converted to float64 anyway by the unit conversion
#         for member in cubelist[1:]:
#             member.coord('time').convert_units(basetimeunit)

def harmonise_time(cubelist):
    "Convert all cubes' time units ('time since...') to that of the first cube in the list."
    basetimeunit=cubelist[0].coord('time').units
    if cubelist[1].coord('time').units != basetimeunit:#Assume the rest aren't right either
        cubelist[0].coord('time').points=cubelist[0].coord('time').points.astype('float64')
        #Because the rest will now be converted to float64 anyway by the unit conversion
        for member in cubelist[1:]:
            member.coord('time').convert_units(basetimeunit)


def harmonise_and_cat(cubelist_in):
    '''Here we need to go through and harmonise attributes that don't match.
    The cubelist is concatenated and returned as a single cube.'''

    #Check for equally lengths of realisation and physics, trim to the shortest
    #(to deal with RCPs having fewer members that historical in CNRM-CM5, for example).
    #This might also be handy if we have to append Short runs to historicals in opatt runs.
    n_realisations=np.empty(len(cubelist_in))
    from iris.exceptions import CoordinateNotFoundError
    #from pdb import set_trace
    #set_trace()

    for i,cube in enumerate(cubelist_in):
        try:
            n_realisations[i]=cube.coord('realization').points.size
        except CoordinateNotFoundError:
            n_realisations[i]=1
    if n_realisations.max()==n_realisations.min():
        cubelist=cubelist_in
    else:
        smallest_cube = cubelist_in[n_realisations.argmin()]
        min_cube_size = n_realisations.min()
        from iris import Constraint
        #print(smallest_cube.coord('realization').points)
        #print(cubelist_in[0].coord('realization').points)
        #for r in cubelist_in[0].coord('realization').points:
        #    print(r in smallest_cube.coord('realization').points)
        trunc_constraint=Constraint(coord_values={'realization':lambda r: r.point in smallest_cube.coord('realization').points})
        from iris.cube import CubeList
        cubelist=CubeList()
        for cube in cubelist_in:
            #realization=cube.coord_dims('realization')[0]
            cubelist.append(cube.extract(trunc_constraint))

    try:
        standard_attr=cubelist[0].attributes
        for member in cubelist[1:]:
            member.attributes=standard_attr
    except AttributeError:
        print('The first cube of this data has no attributes')
    if ( len(cubelist[0].coords('realization'))==1 ):#Check if we have a realization dimension already
        basecoords=cubelist[0].coords(dimensions=cubelist[0].coord_dims('realization'))
        basecoords.extend(cubelist[0].coords(dimensions=[]))#Also gets scalar coordinates
        print(cubelist)
        print(basecoords)
        for member in cubelist[1:]:
            membercoords=member.coords(dimensions=member.coord_dims('realization'))
            membercoords.extend(member.coords(dimensions=[]))#Also gets scalar coordinates
            for coord in [coord for coord in membercoords if coord not in basecoords]:#Remove
                #print 'Removing',coord.long_name
                member.remove_coord(coord)
            for coord in [coord for coord in basecoords if coord not in membercoords]:#Create
                #print 'Adding',coord.long_name
                member.add_aux_coord(coord,cubelist[0].coord_dims(coord.long_name))#member.coord_dims('realization'))
    
    harmonise_time(cubelist)        

    cat=cubelist.concatenate_cube()
    return cat

#def refine_member(enslist,i_memb,membertype,diagstr,forcinglabel):
def refine_member(enslist,i_memb,diagstr,forcinglabel):
    from iris import Constraint
    from iris.util import as_compatible_shape
    from iris.exceptions import ConcatenateError
    #inthismemb=Constraint(cube_func=lambda cube: np.logical_and(cube.attributes[membertype]==i_memb,cube.var_name==diagstr))
    inthismemb=Constraint(cube_func=lambda cube:(cube.attributes['realization']==i_memb//100
                                                 and cube.attributes['physics_version']==i_memb%100
                                                 and cube.var_name==diagstr))
    thismemb=enslist.extract(inthismemb)
    if len(thismemb)>1:#If this member timeseries is split across several files, concatenate them.
        print('Joining '+forcinglabel+'%5d time slices'%i_memb)
#        if i_memb==1303:
#            print(thismemb)
        try:
            member=harmonise_and_cat(thismemb)
        except ConcatenateError:#Likely that some of the cubes have lost their time dimensions. Sort this out!
            from iris.cube import CubeList
            for cube in thismemb:
                for coordname in ['day_of_year','month_number','year']:
                    cube.remove_coord(coordname)#They cause merge problems
            ndims=np.array([len(cube.shape) for cube in thismemb])#The number of dimensions of each cube in the list
            longest=ndims.argmax()#The index of the longest cube
            correctedmemb=CubeList()
            for i,cube in enumerate(thismemb):
                if ndims[i]==ndims[longest]:#just add the cube straight in if it has all the axes
                    correctedmemb.append(cube)
                else:
                    correctedmemb.append(as_compatible_shape(cube,thismemb[longest]))#Rebuild the missing dimension in the cube
            member=harmonise_and_cat(correctedmemb)
            label_dates(member)#Should be safe to put them back now
    else:
        member=thismemb[0]
        print('Refined '+forcinglabel+'%5d time slices'%i_memb)
    repair_coord_system(member)#Repair issue with Iris 1.7 regridding
    return member


def harmonise_calendars(obs,all_ens):
    'Convert obs onto the same calendar as the model (usually the better way around for CPU and memory!)'
    from iris.cube import CubeList
    from cf_units import Unit
    #from scipy.interpolate import interp1d
    obslist=CubeList()
    thisyear=obs.coord('year').points.min()
    finalyear=obs.coord('year').points.max()
    dropdays=np.array([76,132,188,243,298])#17 Mar, 12 May, 7 Jul, 31 Aug, 25 Oct
    while thisyear<=finalyear:
        obsthisyear=obs[obs.coord('year').points==thisyear]
        alltimeyear=all_ens.coord('time')[all_ens.coord('year').points==thisyear]
        alltime_orig=alltimeyear.copy()#So we have the original time units to hand
        #Now convert both time units into years since the beginning of this year
        timeunitname='day since %4d-01-01 00:00:00.0000000 UTC'%thisyear#Couldn't we just do this with a day_of_year AuxCoord?
        if alltimeyear.units.name != timeunitname:
            alltimeyear.convert_units(Unit(timeunitname,calendar=alltimeyear.units.calendar))#Does this crash if these are already the units?
        if obsthisyear.coord('time').units.name != timeunitname:
            obsthisyear.coord('time').convert_units(Unit(timeunitname,calendar=obsthisyear.coord('time').units.calendar))
        #time_interp=interp1d(obsthisyear.coord('time').points,obsthisyear.data,axis=0)
        #newobsdata=time_interp(alltimeyear.points)#Interpolate onto model time axis
        #interp1d can't deal with masks. Possibly a flawed idea anyway.
        #Instead, drop 29 Feb (if leap), 17 Mar, 12 May, 7 Jul, 31 Aug, 25 Oct so 90 days in each season.
        if len(obsthisyear.coord('time').points)==366:
            dropthisyear=np.hstack([59,dropdays])
        else:
            dropthisyear=dropdays-1#Zero indexed, remember!
        newobsdata=np.ma.array( np.delete(obsthisyear.data.data, dropthisyear, 0), mask=np.delete(obsthisyear.data.mask, dropthisyear, 0) )
        #Now rebuilt obsthisyear with new time axis and interpolated data
        obsthisyear.remove_coord('time')
        obsthisyear=obsthisyear[0:len(alltime_orig.points)]#Truncate cube to right length
        obsthisyear.data=newobsdata#Replace data with interpolated data
        obsthisyear.add_dim_coord(alltime_orig,0)#Replace obs time axis with model time axis
        
        obslist.append(obsthisyear)
        thisyear+=1
    return obslist.concatenate_cube()


def mask_ens(cubelist,mask_in):
    '''Masks an ensemble cubelist. For an obs cube use mask_in=obs.data.mask. Should also preserve any masking already present'''
    if isinstance(cubelist[0].data,np.ma.core.MaskedArray):#Type check to see if the data is already masked.
        print('Combining obs and model masks')
        for member in cubelist:#Combine the masks - if it's masked in either of the old arrays, it's still masked now.
            member.data.mask=np.logical_or(mask_in,member.data.mask)
    else:#It's not masked already, so make it masked.
        print('Applying obs mask to model')
        for member in cubelist:
            member.data=np.ma.array(member.data,mask=mask_in)


def mask_cube(cube,mask_in,maskdims=3):
    '''Masks an ensemble cube. For an obs cube use mask_in=obs.data.mask. Should also preserve any masking already present'''
    tilesize=np.array(cube.shape)
    tilesize[-maskdims:]=1#Set the last 3 dimensions to not multiply up (mask is already the right length and shape)
    tiled_mask=np.tile(mask_in, tilesize)#Repeat mask across all members
    #tiled_mask=np.tile(mask_in, [cube.coord('realization').points.size,1,1,1])#Repeat mask across all members
    if isinstance(cube.data,np.ma.core.MaskedArray):#Type check to see if the data is already masked.
        print('Combining obs and model masks')
        #Combine the masks - if it's masked in either of the old arrays, it's still masked now.
        cube.data.mask=np.logical_or(tiled_mask,cube.data.mask)
    else:#It's not masked already, so make it masked.
        print('Applying obs mask to model')
        print(cube.data.shape, tiled_mask.shape)
        cube.data=np.ma.array(cube.data,mask=tiled_mask)

def get_landsea_mask(lfpath,egcube,masksea=True):
    '''From a landfraction file, makes a mask where not land==True
     on the same grid as an example cube in a given Giorgi region.'''
    from iris import load_cube
    landfraction_raw=load_cube(lfpath)
#    areaconstraint=iris.Constraint(coord_values={'latitude': lambda lat: coord[1]<=lat<=coord[3],
#                                                 'longitude': lambda lon:coord[0]<=lon<=coord[2]})
    repair_coord_system(landfraction_raw)#Repair issue with Iris 1.7 regridding
    lf=landfraction_raw.regrid(egcube,ia.AreaWeighted(mdtol=0.333))#.extract(areaconstraint)#Just use constraints already on egcube
    if masksea:
        mask=(lf.data<80)#Land fraction of less than 80% is not land (aka notlandmask)
    else:
        mask=(lf.data>20)#Land fraction of greater than 20% is not sea (aka notseamask)
    return mask


def coord_or_attrib_value(cube,name):
    'Return the value of a given variable which could be stored in either a coordinate or an attribute.'
    if len(cube.coords(name))>0:
        value=cube.coord(name).points
    else:
        value=cube.attributes[name]
    return value


def seasonal_ensemble(cube_ens,startmonth,endmonth,obs_seas_timeaxis):
    '''Makes seasonal average of each ensemble member, for the derived start and end month of the season.
    Also harmonises the time axis with that input from the obs for the season.'''
    from iris.exceptions import CoordinateCollapseError
    try:
        ens_seas=filter_months(cube_ens,startmonth,endmonth).aggregated_by('year',ia.MEAN)
        ens_seas.coord('time').points=obs_seas_timeaxis.points
        ens_seas.coord('time').bounds=obs_seas_timeaxis.bounds
        ens_seas.coord('time').units=obs_seas_timeaxis.units
    except CoordinateCollapseError:#This occurs when there's only 1 months and 1 year, e.g. seasonal runs
        ens_seas=filter_months(cube_ens,startmonth,endmonth)#.collapsed('time',ia.MEAN)
        #Will we still need the time axis though, just with 1 point?
    return ens_seas


def delayed_loop(delay_options,delayed_function,*args,**kwargs):
    '''Useful for ftp puts and MASS gets, where we might have to wait for the server to come back up.'''
    for delay in delay_options:
        try:
            output=delayed_function(*args,**kwargs)
            break#if successful, quit this innermost loop - no further delays required
        except IOError:
            delaymessage='Function failed, waiting %1d hour'%delay
            if delay>1:
                delaymessage+='s'
            if delay==-1:
                print("Max delay expired, giving up on running the function.")
                return
            else:
                print(delaymessage)
                sleep(delay*3600)#seconds
    else:#Just in case the -1 doesn't trigger exit, we have an else at the end of the for loop
        return#If no break occurs (file retrieval was never successful), don't continue the loop nest
    return output#This should be the only place the data gets out


'Here follows code for optimal detection.'


def icom2inv_trans_lin(icom):
    '''A port of the section of opt_detn.pro which converts icom into an inverted matrix.
    The digit selection is now done with integer division instead of string operations,
    the total digit count by base-10 logarithm, and the j loop is vectorised.'''
    icomlen=len(icom)#nez in opt_detn.pro
    nforcings=np.floor(np.log10(icom))+1#nedd in opt_detn.pro - This is the number of digits of each element
    lin_matrix=np.matrix(np.zeros([icomlen,icomlen],dtype='int32'))
    for i,icom_element in enumerate(icom):
        jforcings=np.arange(nforcings[i],0,-1)#Work our way through the digits of icom_element from left to right
        k=np.int32(icom_element%10**jforcings#Remove digits ahead of the one of interest
                            //10**(jforcings-1))#and behind to give the value of each digit of forcing in an icom element
        lin_matrix[i,k-1]=1#IDL references in opposite order to Python. This means we don't need the later transpose.
    inv_trans_lin=lin_matrix.I#Invert the matrix
    return inv_trans_lin


def icom2inv_trans_lin_bool(icom):
    'Alternative boolean version of icom2inv_trans_lin'
    icomlen=len(icom)#nez in opt_detn.pro
    flabels=np.tile(np.arange(icomlen)+1, [icomlen,1])#List of potential forcing digits, repeated over icom length
    icomtile=np.tile(icom,[icomlen,1]).transpose()#icom made square for comparison operation
    jforcing=icomlen-1
    lin_matrix_v=(icomtile//10**jforcing==flabels)#Initialises matrix by checking presence left-most digit of longest possible icom listing
    while jforcing>0:
        lin_matrix_v=np.logical_or(lin_matrix_v, icomtile%10**jforcing              #Now check the next digit down from previous down
                                                        //10**(jforcing-1)==flabels)#and combine it (OR) with tally for previous digits
        jforcing-=1#Work our way down in magnitude through the digits til we've checked and tallied them all in the matrix
    inv_trans_lin=np.matrix(lin_matrix_v, dtype='int32').I#Convert and invert the matrix
    return inv_trans_lin


def nparr2Rmatrix(x,transpose=False):
    '''Adapted from code by someone called Tim on
    http://stackoverflow.com/questions/14463600/converting-numpy-array-to-rpy2-matrix-forecast-package-xreg-parameter'''
    import rpy2.robjects as R

    if np.logical_or(transpose==True, x.shape[0]==1):
        #print 'Transposing input array'
        xvec = R.FloatVector(x.transpose().reshape((x.size,1)))
        nc, nr = x.shape
    else:
        xvec = R.FloatVector(x.reshape(x.size,1))
        nr, nc = x.shape
    xr = R.r.matrix(xvec, nrow=nr, ncol=nc, byrow=True)#, dimnames=dimnames)
    return xr


def plot_histogram(dataset,title='Histogram',xlabel='Data',ylabel='Frequency'):
    'Really generic histogram wrapper'
    from matplotlib.pyplot import figure,hist,show
    fig=figure(num=title)
    subpl=fig.add_subplot(1,1,1)
    subpl.set_xlabel(xlabel)
    subpl.set_ylabel(ylabel)
    hist(dataset)
    show(block=False)
       

#def map_stats(stats,statname,show_graph=True,filename='map_stats.png',symmetry=True):
#    from iris.plot import pcolormesh
#    from matplotlib.pyplot import figure,gca,colorbar,show,close
#    from matplotlib.cm import get_cmap
#    fig=figure(num='Map of '+statname)
#    subpl=fig.add_subplot(1,1,1)
#    if (str(stats.units)[:2]=='mm') or (str(stats.units)[:6]=='kg m-2'):#Or actually compare units?
#        colourmap=get_cmap('BrBG')
#    else:
#        colourmap=get_cmap('bwr')
#    #colourmap.set_bad('#c0c0c0')#Grey for missing data, to distinguish from white for no change.
#    if symmetry:
#        #Now make a symmetrical scale to apply this colourmap to
#        #absmax=np.ma.absolute(stats.data).max()
#        absmax=3.*stats.data.std()
#        pcolormesh(stats,vmin=-absmax, vmax=absmax, cmap=colourmap)
#    else:#just let pcolormesh do what it likes
#        pcolormesh(stats,cmap=colourmap)
#    gca().coastlines()
#    colorbar(orientation='vertical',label=str(stats.standard_name)+'('+str(stats.units)+')')
#    if show_graph:
#        show(block=False)
#    else:
#        fig.savefig(filename)
#        close(fig)


def map_subplot(stats,symmetry,vrange=None):
    from iris.plot import pcolormesh,pcolor
    from matplotlib.pyplot import gca#,colorbar
    from matplotlib.cm import get_cmap
    if (str(stats.units)[:2]=='mm') or (str(stats.units)[:6]=='kg m-2'):# or (str(stats.units)[3:9]=='kg m-2'):#Or actually compare units?
        colourmap=get_cmap('BrBG')
    else:
        #colourmap=get_cmap('bwr')
        colourmap = get_cmap('coolwarm')
    #colourmap.set_bad('#c0c0c0')#Grey for missing data, to distinguish from white for no change.
    #colourmap.set_bad('#000000',alpha=0)#Black
    if type(vrange) is list:
        print('Using specified range')
        pcolormesh(stats, vmin=vrange[0], vmax=vrange[1], cmap=colourmap)
    elif symmetry:
        from pdb import set_trace
        set_trace()
        #Now make a symmetrical scale to apply this colourmap to
        #absmax=np.ma.absolute(stats.data).max()
        absmax=3.*stats.data.std()
        pcolormesh(stats,vmin=-absmax, vmax=absmax, cmap=colourmap)#, hatch='/')
    else:#just let pcolormesh do what it likes
        pcolormesh(stats,cmap=colourmap)
        pcolor(stats,alpha=0, hatch='/')
    gca().coastlines()


def map_stats(stats,statname,subnames=None,show_graph=True,filename='map_stats.png',symmetry=True,vrange=None):
    from matplotlib.pyplot import figure,show,close,colorbar
    fig=figure(num='Map of '+statname)
    statdims=stats.shape
    if len(statdims)==2:
        subpl=fig.add_subplot(1,1,1)
        map_subplot(stats,symmetry,vrange=vrange)
        colorbar(orientation='vertical', label=str(stats.long_name) + '(' + str(stats.units) + ')')
    elif len(statdims)==3:
        subplcols=np.ceil(np.sqrt(statdims[0]))
        subplrows=np.ceil(statdims[0]/subplcols)
        for i in range(statdims[0]):
            if type(subnames)==list:
                subtitle=subnames[i]
            else:
                subtitle=''
            subpl=fig.add_subplot(subplrows,subplcols,i+1,title=subtitle)
            map_subplot(stats[i],symmetry,vrange=vrange)
        colorbar(orientation='vertical',label=str(stats[0].long_name)+'('+str(stats[0].units)+')')
    if show_graph:
        show(block=False)
    else:
        fig.savefig(filename)
        close(fig)


    
def plot_distn(Taxis,distn,varn,label,obs=None,comparator='gt',hatch=None,memberlines=False,**kwargs):
    'Calculate and plot a Gaussian around the given distribution.'
    #import matplotlib.pyplot as mp
    from matplotlib.pyplot import plot, fill_between
    if varn=='pr':
        import scipy.stats as ss
        fit_alpha,fit_loc,fit_beta=ss.gamma.fit(distn,loc=0)
        probden=ss.gamma.pdf(Taxis,fit_alpha,fit_loc,fit_beta)
        #shape,rate=r_fit_gamma(distn)
        #probden=r_gamma_distn(Taxis,shape,rate)
        if obs==None:
            peak=probden.mean()*1.5
        else:
            peak=ss.gamma.pdf(obs,fit_alpha,fit_loc,fit_beta)
    elif varn=='Risk Ratio':
        from warnings import catch_warnings,simplefilter
        with catch_warnings():
            simplefilter("ignore")#UserWarning: invalid value encountered in log2
            #print('log Gaussian')
            log2rr=np.log2(distn)
            probden=gaussian_distn(np.log2(Taxis),log2rr.mean(),log2rr.std())
            peak=None
    else:#elif varn=='tas':
        # probden=gev.pdf(Taxis,*gev.fit(distn)) # GEV fit could be done like this
        Tmean=distn.mean()
        Tstdev=distn.std()
        probden=gaussian_distn(Taxis,Tmean,Tstdev)
        if obs==None:
            peak=probden.mean()*1.5
        else:
            peak=gaussian_distn(obs,Tmean,Tstdev)
    #peak=probden.mean()*1.5
    graph=plot(Taxis,probden,label=label,**kwargs)#Plot Gaussian
    if varn=='Risk Ratio':
        graph.append(fill_between(Taxis,probden,0,hatch=hatch,**kwargs))
    elif obs!=None:#Shade beyond threshold
        if comparator=='gt' or comparator=='ge':
            Taxis_ge_obs=(Taxis>=obs)
            graph.append(fill_between(np.append(obs,Taxis[Taxis_ge_obs]),
                                      np.append(peak,probden[Taxis_ge_obs]),0,hatch=hatch,**kwargs))#
        elif comparator=='lt' or comparator=='le':
            Taxis_le_obs=(Taxis<=obs)
            graph.append(fill_between(np.append(Taxis[Taxis_le_obs],obs),
                                      np.append(probden[Taxis_le_obs],peak),0,hatch=hatch,**kwargs))#
    if len(distn)<=10 or memberlines==True:
        for membval in distn:
            #Plot each actual value as a short vertical line in the distribution (remove this loop if too many).
            graph.append(plot(np.array([membval,membval]),np.array([0,0.2/np.abs(distn).mean()]),**kwargs))
    return graph,peak



def TREND(data,axis,ordinates):
    'Custom aggregator for gridpoint trends, using linear algebraic least-squares.'
    'Calculation as described in documentation for np.linalg.lstsq().'
    'Further ideas lifted from scipy.signal.detrend().'
    size=data.shape
    if axis<0:#Make axis number positive
        axis+=len(size)
    #Now, transpose data so the axis we want is the first dimension
    newshape=np.r_[axis,0:axis,axis+1:len(size)]
    newsize=np.hstack([ size[axis], size[0:axis], size[axis+1:len(size)] ]).astype(np.int32)
    dt=data.transpose(newshape).reshape([newsize[0],newsize[1:].prod()])#data must be a 2d matrix
    #A=np.vstack([ordinates,np.ones(len(ordinates))]).T #Ones column accounts for an offset
    #trend=np.linalg.lstsq(A,dt)[0][0] #Trim off all the other info that lstsq gives us, eg. offset
    #The above is a problem because it doesn't deal with the masked data correctly. Do it gridboxwise instead.
    trendarr=np.ma.empty(newsize[1:].prod())
    i=0
    while i<trendarr.size:
        not_missing=np.logical_not(dt[:,i].mask)#Which years have data
        n_not_missing=not_missing.sum()
        if n_not_missing>1:#Need at least 2 points to make a trend!
            A=np.vstack([ordinates[not_missing],np.ones(n_not_missing)]).T
            trendarr[i]=np.linalg.lstsq( A, dt[:,i].data[not_missing] )[0][0]
        i+=1
    trendarr.mask=(dt.mask.sum(0)>=ordinates.size/2.)#Additional filter for how many missing years we'll allow. 50% missing tolerance here
    return trendarr.reshape(newsize[1:])
    #return np.ma.array( trendarr.reshape(newsize[1:]), mask=(data.mask.sum(axis)==0) )#Mask if any time point was masked - do better!
    #Put trend matrix back into the original array shape and return


def coord_diffs(coord1, coord2):
    """
    Ruth Comer's function. Print which aspects of two coords do not match.
    """
    for attr in ['standard_name', 'long_name', 'var_name', 'units',
                 'coord_system', 'circular']:
        if getattr(coord1, attr, None) != getattr(coord2, attr, None):
            print('{} difference: {} vs {}'.format(
                attr, repr(getattr(coord1, attr, None)),
                repr(getattr(coord2, attr, None))))

    common_keys = set(coord1.attributes).intersection(coord2.attributes)
    for key in set(coord1.attributes) - common_keys:
        print('attribute only on coord1: {}'.format(key))
    for key in set(coord2.attributes) - common_keys:
        print('attribute only on coord2: {}'.format(key))
    for key in common_keys:
        if np.any(coord1.attributes[key] != coord2.attributes[key]):
            print('attributes differ: {}'.format(key))

    for attr in ['points', 'bounds']:
        array1 = getattr(coord1, attr)
        array2 = getattr(coord2, attr)
        if (array1 is None) != (array2 is None):
            has_array, = [coord for (coord, array) in
                          zip(['coord1', 'coord2'], [array1, array2]) if
                          array is not None]
            print('{} exist only on {}'.format(attr, has_array))
        else:
            if np.any(array1 != array2):
                print('{} differ'.format(attr))
            if array1 is not None and array1.dtype != array2.dtype:
                print('{} dypes differ: {} vs {}'.format(
                    attr, array1.dtype, array2.dtype))
