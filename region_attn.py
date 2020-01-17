#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# No longer #!/usr/local/sci/bin/python2.7
'''
Created on May 21, 2014

@author: Fraser Lott, Met Office Hadley Centre. Crown Copyright.
'''
import numpy as np
import iris
import iris.analysis as ia
import fl_detn_tools as fl

# from scipy.stats import genextreme as gev
# import matplotlib.pyplot as plt
# import matplotlib.pyplot as plt

clim_period = np.array([1961,
                        1990])  # This should usually be set to the climatology period of the obs anomalies, or alternatively apply this to the obs too.


def distributed_load(*args, **kwargs):
    'Informatics Lab parallelised loading via Dask, for use with Iris 2.'
    from dask import delayed
    from dask.bag import from_delayed
    from glob import glob

    file_list = glob(args[0])

    # print(file_list)

    @delayed
    def load_s(fname):
        return iris.load_raw(fname, *args[1:], **kwargs)

    cubes = from_delayed([load_s(f) for f in file_list])
    results = iris.cube.CubeList(cubes.compute())
    return results.merge()


def region_constraint(coord):
    'This version uses year pre-labelled by fl.label_dates(). Comments show alternative.'
    # coord=fl.giorgi_coord(region)
    # Allow for generic regions by moving search for Giorgi coord outside.

    #    coord=fl.trieste_regions_coord(region)
    #    import datetime as dt
    #    timeunit=cube.coord('time').units
    #    starttime=timeunit.date2num(dt.datetime(startyear,1,1,0,0,0))
    #    endtime=timeunit.date2num(dt.datetime(endyear,12,31,23,59,59))
    return iris.Constraint(coord_values={'latitude': lambda lat: coord[1] <= lat <= coord[3],
                                         'longitude': lambda lon: coord[0] <= lon <= coord[2]})


#                                         'year': lambda year: startyear<=year<=endyear})
#                                         'time': lambda date: starttime<=date<=endtime})
# With year labels, this was much easier - otherwise we need to make datetime objects for comparison

def mask_obs(obs, varn, varo, freq, notlandmask):
    if varn == 'tas':
        # Combine the obs mask with the land mask. This will then be propagated to the models too.
        obs.data = obs.data.copy()#Should make mask writable in numpy>=1.16
        obs.data.mask = np.logical_or(obs.data.mask, notlandmask)
    elif varn == 'pr':
        print('Redefining units')
        if varo == 'Ann' or freq == 'day':  # GHCND
            obs.units = 'kg m-2 day-1'
        else:
            obs.units = 'kg m-2 month-1'
        bigmask = np.tile(notlandmask.reshape(1, obs.shape[1], obs.shape[2]), (obs.shape[0], 1, 1))
        obs.data = np.ma.array(obs.data, mask=bigmask)
    elif varn == 'mrsos':
        print('Redefining units')
        # obs.units.convert('m')
        obs.units = 'kg m-2'  # Thanks to the wonders of SI units


def calc_ensemble_anomaly(ensemble, clim_period):
    'Turns actuals into anomalies for each member of an ensemble cubelist.'
    ens_anom = iris.cube.CubeList()
    for member in ensemble:
        ens_anom.append(fl.anomalise(member, clim_period[0], clim_period[1]))
    return ens_anom


def obs_and_landmask_anomalise_model(model_ens, varn, obs, clim=None):
    if model_ens.units != obs.units:
        # if model_ens[0].units != obs.units:
        # for i,_ in enumerate(model_ens):
        #    model_ens[i].convert_units(obs.units)
        model_ens.convert_units(obs.units)
    if (varn == 'tas') | (varn == 'psl'):
        # model_anom=calc_ensemble_anomaly(model_ens,clim_period)
        model_anom = fl.anomalise(model_ens, clim_period[0], clim_period[1], clim=clim)
        #print('No anomalies here')
    else:  # elif varn=='pr':
        model_anom = model_ens
    # fl.mask_ens(model_anom,obs.data.mask)
    fl.mask_cube(model_anom, obs.data.mask)
    if (varn == 'mrso'):  # Catch issues with stray negative points which didn't get masked
        # for i,member in enumerate(model_ens):
        #    model_anom[i].data.mask[member.data.data<0]=True
        model_anom.data.mask[model_anom.data.data < 0] = True
    return model_anom  # ,Changes to obs masking should pass though automatically


def get_savefile_names(forcinglabel, startyear, endyear, region, modelroot, modelname, diagstr, freq, membertype,
                       uniquemembs):
    #    regiontitle_dict={"<class 'ascend.shape.Shape'>":lambda region: region.attributes['name'],
    #                      "<type 'list'>":lambda region: region[4],
    #                      "<type 'bool'>":lambda region: ''}#For region=False (i.e. global mean)
    #    regiontitle=regiontitle_dict[str(type(region))](region)#Choose region title from the above dictionary

    # from ascend.shape import Shape #Temporarily comment out this and the line below with Ascend problems
    regiontitle_dict = {  # Shape:lambda region: region.attributes['name'],
        list: lambda region: region[4],
        bool: lambda region: ''}  # For region=False (i.e. global mean)
    regiontitle = regiontitle_dict[type(region)](region)  # Choose region title from the above dictionary
    saveroot = modelroot + '/saves/' + diagstr + '_' + freq + '_' + modelname + '_' + forcinglabel + '_r'
    savetail = '_%4d01_%4d12_' % (startyear, endyear) + regiontitle + '.nc'
    if membertype:  # is present, make a list of filenames (assumes uniquemembs also specified)
        if membertype == 'physics_version':
            saveroot += '1i1p'
        elif membertype == 'realization':
            savetail = 'i1p1' + savetail
        savenames = []
        for i_memb in uniquemembs:
            savenames.append(saveroot + str(i_memb) + savetail)
    else:  # Just do a wildcarded filename
        savenames = saveroot + '*' + savetail
    return savenames


def load_region_ensemble(forcinglabel, region, startyear, endyear, obs, modelroot, modelname, diagstr, freq,
                         save=False,startendmonths=None):
    'Loads entire model ensemble for specified time period and region, regridded to observational grid.'
    # iris.FUTURE.cell_datetime_objects=True#so that this will behave the same in Iris 1 or 2.
    # Doesn't behave the same between Python 2 and 3 though!
    # from sys import version#Get Python version for checking.
    # This assumes we have Iris 2 with Python 3, or Iris 1 with Python 2.

    #    if modelname=='HadGEM3-A-N216':
    #        membertype='physics_version'#Set for perturbed physics ensembles
    #        memberwild='_r1i1p*'#'_r*i1p*'#
    #    else:
    #        membertype='realization'#Set for initial condition ensembles
    #        memberwild='_r*i1p1*'
    memberwild = '_r*i1p*'
    if modelroot == '/project/detn/' or modelroot == '/s3/informatics-eupheme/':  # +modeln+'/historical/'+varn+'/'+freq+'/'
        modelpath = modelroot + modelname + '/' + forcinglabel + '/' + diagstr + '/' + freq + '/' + diagstr + '_' + freq + '_' + modelname + '_' + forcinglabel + memberwild + '.nc'
    else:
        modelpath = modelroot + '/' + modelname + '/' + diagstr + '_' + freq + '_' + modelname + '_' + forcinglabel + memberwild + '.nc'
    if type(startendmonths)==list:
        #TODO: alter dateconstraint to account for seasons crossing year boundary.
        # For now we leave in all the months in that situation and deal with it later
        if startendmonths[0]>startendmonths[1]:
            dateconstraint = iris.Constraint(coord_values={'year': lambda year: startyear <= year <= endyear+1})
        else:
            dateconstraint = iris.Constraint(coord_values={'year': lambda year: startyear <= year <= endyear,
                                                           'month_number': lambda mon: startendmonths[0]<=mon<=startendmonths[1]})
    elif forcinglabel != 'piControl':
        dateconstraint = iris.Constraint(coord_values={'year': lambda year: startyear <= year <= endyear})
    #    enslist=iris.load(modelpath,areaconstraint&dateconstraint,callback=fl.datelabel_callback)#Read in all model ensemble members matching the wildcards in modelpath
    # else:
    #    enslist=iris.load(modelpath,areaconstraint,callback=fl.datelabel_callback)#Read in all model ensemble members matching the wildcards in modelpath
    from warnings import catch_warnings, simplefilter
    with catch_warnings():
        simplefilter("ignore")
        # UserWarning: Missing CF-netCDF measure variable u'areacella', referenced by netCDF variable u'tas'
        if iris.__version__.startswith('2'):  # then we can use Dask for parallel loading
            enslist = distributed_load(modelpath, dateconstraint, callback=fl.datelabel_callback)
        else:  # we'll load it all serially
            enslist = iris.load(modelpath, dateconstraint,
                                callback=fl.datelabel_callback)  # Read in all model ensemble members matching the wildcards in modelpath
    # Area constraint is not all that necessary, as this is propagated from the obs when we regrid
    # However, it might be needed for processing dailies to keep memory levels reasonable.

    membnum = np.empty(len(enslist), dtype='int32')
    for i, member in enumerate(enslist):
        # membnum[i]=member.attributes[membertype]#Read this in case we read the members in the wrong order
        r = fl.coord_or_attrib_value(member, 'realization')
        p = fl.coord_or_attrib_value(member, 'physics_version')
        membnum[i] = r * 100 + p  # Maybe not pythonic, but works with np.unique
        # membnum[i]=member.attributes['realization']*100+member.attributes['physics_version']

    uniquemembs = np.unique(membnum)
    # Get the length of the time axis
    member_constrained = apply_region(region, fl.refine_member(enslist, uniquemembs[0], diagstr, forcinglabel), diagstr)
    #if iris.__version__.startswith('2.0'):  # where datetime has comparators and Iris uses datetime
    #    endtime = max(member_constrained.coord('time').cells())
    #    constr_time1st = iris.Constraint(coord_values={'time': lambda time: time <= endtime})
    #from pdb import set_trace
    #set_trace()
    if iris.__version__.startswith('1'):  # python2.7, Iris 1. datetimes aren't comparable and Iris uses Julian days
        endtime = member_constrained.coord('time').points.max()
        constr_time1st = iris.Constraint(coord_values={'time': lambda time: time <= endtime})
    else:
        time_unit = member_constrained.coord('time').units
        endtime = time_unit.num2date(member_constrained.coord('time').points.max())
        constr_time1st = iris.Constraint(coord_values={'time': lambda time: time.point <= endtime})
    ##else:  # Iris 2.1 using cftime datetimes which are incomparable again. Need to key to cells() to points for max, then get a PartialDateTime
    ##    endtime_cell = max(member_constrained.coord('time').cells(), key=lambda cell: cell.point)
        # from iris.time import PartialDateTime
        # endtime=PartialDateTime(year=endtime_cell.point.year, month=endtime_cell.point.month, day=endtime_cell.point.day)
        # yearconstr=iris.Constraint(coord_values={'year':lambda year: year<endtime_cell.point.year})
    ##    constr_time1st = iris.Constraint(time=lambda cell: np.logical_or(cell.point.year < endtime_cell.point.year,
    ##                                                                     np.logical_and(
    ##                                                                         cell.point.year == endtime_cell.point.year,
    ##                                                                         cell.point.dayofyr <= endtime_cell.point.dayofyr)))

    # constr_time1st=iris.Constraint(coord_values={'time': lambda time: time<=endtime})#starttime<=time
    def subsel_regrid(i, i_memb):  # Nested function for parallel load on Iris 2
        member = fl.refine_member(enslist, i_memb, diagstr, forcinglabel)
        member_constrained = apply_region(region, member, diagstr)
        fl.prep_for_merge(member_constrained)
        if i > 0:
            # Check the time axis isn't longer than the first member, trim if necessary
            member_constrained = member_constrained.extract(constr_time1st)
        return member_constrained.regrid(obs, ia.AreaWeighted(mdtol=0.333))  # Regrid onto obs grid

    from psutil import virtual_memory
    if iris.__version__.startswith('2') and virtual_memory()._asdict()['total']>8e9:  # Then we can use delayed distributed load
        from dask import delayed
        delayed_regrid = delayed(subsel_regrid)
        list_of_cubes = [delayed_regrid(i, i_memb) for i, i_memb in enumerate(uniquemembs)]
        regionens = delayed(iris.cube.CubeList)(list_of_cubes).compute()
    else:  # Iris 1
        list_of_cubes = [subsel_regrid(i, i_memb) for i, i_memb in enumerate(uniquemembs)]
        regionens = iris.cube.CubeList(list_of_cubes)

    if save:  # This will probably need updating if we use full rip
        savename = get_savefile_names(forcinglabel, startyear, endyear, region, modelroot, modelname, diagstr, freq,
                                      membertype, uniquemembs)
        for i in range(len(regionens)):
            iris.save(regionens[i], savename[i])
        print('Exiting...')
        from sys import exit
        exit()  # Only temporary! This is a really bad way to kill the program
    else:
        #print(regionens)
        #from pdb import set_trace
        #set_trace()
        from iris.exceptions import MergeError
        try:
            return regionens.merge_cube()
        except MergeError:
            fl.harmonise_time(regionens)#In case some cubes have a different time base (can happen with GISS-E2-H)
            return regionens.merge_cube()


def load_and_regrid_cube(forcinglabel, region, startyear, endyear, obs, modelroot, modelname, diagstr, freq):
    'Reproduces load_region_ensemble() for when the ensemble data is already saved in a global cube.'
    modelpath = modelroot + '/' + modelname + '/' + diagstr + '_' + freq + '_' + modelname + '_' + forcinglabel + '.nc'
    enscube = iris.load_cube(modelpath)
    constrained = apply_region(region, enscube, diagstr)
    return constrained.regrid(obs, ia.AreaWeighted(mdtol=0.333))  # Regrid onto obs grid


def constrain_shape(region, obs_raw, varname_o):
    return region.constrain_cube(obs_raw)


def constrain_box(region, obs_raw, varname_o):
    "Otherwise assume region is specified by a list of coordinates and use it as a constraint."
    areaconstraint = region_constraint(region)  # Choose region and time period
    return obs_raw.extract(areaconstraint)


def constrain_box_with_extent(region, obs_raw, varname_o):
    "Uses intersection instead of constraint to deal with circular coords (particularly longitude)."
    obs_raw.coord('longitude').circular = True  # Even if GPCC says otherwise
    constrained = obs_raw.intersection(longitude=(region[0], region[2]), latitude=(region[1], region[3]),
                                       ignore_bounds=True)
    for coordname in ['latitude', 'longitude']:
        if not constrained.coord(coordname).is_contiguous():
            constrained.coord(coordname).bounds = None
            constrained.coord(coordname).guess_bounds()
    return constrained


def no_constraint(region, obs_raw, varname_o):
    return obs_raw


def apply_region(region, obs_raw, varname_o):
    # from ascend.shape import Shape #Temporarily comment out this and the line below with Ascend problems
    applydict = {  # Shape:constrain_shape,#e.g. from fl.get_srex_shape()
        list: constrain_box_with_extent,  # Giorgi region corners
        bool: no_constraint}  # For region=False (i.e. global mean)
    region_obs = applydict[type(region)](region, obs_raw, varname_o)
    return region_obs


def tamsat_callback(cube, field, fname):
    cube.coord('time').convert_units('days since 1983-01-01')  # or months
    for metatitle in ['file01', 'file02', 'file03']:
        metadata = cube.attributes.pop(metatitle)
        cube.add_aux_coord(iris.coords.AuxCoord(metadata, long_name=metatitle),
                           data_dims=0)  # Associate with time dimension
    _ = cube.attributes.pop('history')
    fl.label_dates(cube)


def ghcnd_callback(cube, field, fname):
    'Build a better time unit from the yyyymmdd (with incorrect mm!) given in GHCND NetCDF files'
    from datetime import datetime
    from time import strptime
    from cf_units import Unit
    timeunit = Unit('days since 1950-01-01 00:00:00', calendar='gregorian')
    year = cube.coord('time').points // 1e4
    if cube.var_name != 'Ann':
        month = strptime(cube.var_name,
                         '%b').tm_mon  # The month is wrong in the time dimension, but is right in the variable name (!!)
        cube.var_name = 'pr'
    else:
        month = 1
    # month=cube.coord('time').points%1e4//1e2
    # day=cube.coord('time').points%1e2
    cube.remove_coord('time')

    timedata = np.empty(len(year), dtype='float64')
    for i, _ in enumerate(timedata):  # Can't seem to vectorise datetime objects
        # timedata[i]=timeunit.date2num(datetime(int(year[i]),int(month[i]),int(day[i]),0,0,0))#This date is wrong in the file!!
        timedata[i] = timeunit.date2num(datetime(int(year[i]), month, 1, 0, 0, 0))
    newtime = iris.coords.DimCoord(timedata, standard_name='time', units=timeunit)
    cube.add_dim_coord(newtime, 0)

    # fl.label_dates(cube)
    from iris.coord_categorisation import add_year
    add_year(cube, 'time', name='year')  # Don't add month_number yet, else it gets merged later as a new dimension


def ghcnd_cat(obs_list):
    "GHCND time slices are all out of order, so won't concatenate. The following rearranges the order."
    obs_full = iris.cube.CubeList()
    n_tslices = obs_list[0].shape[0]
    i = 0
    while i < n_tslices:
        for cube in obs_list:
            obs_full.append(cube[i])
        i += 1
    obs_cube = obs_full.merge_cube()
    from iris.coord_categorisation import add_month_number
    add_month_number(obs_cube, 'time', name='month_number')
    return obs_cube


def load_obs(startyear, endyear, region, obsset, varname_o, daily=False):  # ,obsname):
    '''Loads an observations cube from the NetCDF file.'''
    obsname, varname_o, _, _ = fl.set_variables(obsset, varname_o, daily=daily)
    dateconstraint = iris.Constraint(coord_values={'year': lambda year: startyear <= year <= endyear})
    if obsset == 'GHCNDEX_Rx1day':  # (varname_o=='Ann'):
        nameconstraint = iris.Constraint(
            cube_func=lambda cube: cube.var_name != varname_o)  # Get each month but not annual for GHCND
    else:
        nameconstraint = iris.Constraint(cube_func=lambda cube: cube.var_name == varname_o)
    if (obsset == 'GHCNDEX_Rx1day') | (obsset == 'TAMSAT') | (
            obsset == 'aphrodite'):  # (varname_o=='rfe') | (varname_o=='Ann'):#Tamsat rainfall estimates and GHCND files need concatenating
        callback_dict = {'TAMSAT': tamsat_callback, 'GHCNDEX_Rx1day': ghcnd_callback,
                         'aphrodite': fl.datelabel_callback}
        if iris.__version__.startswith('2'):
            obs_list = distributed_load(obsname, nameconstraint & dateconstraint, callback=callback_dict[obsset])
        else:
            obs_list = iris.load(obsname, nameconstraint & dateconstraint, callback=callback_dict[obsset])
        if obsset == 'aphrodite':
            fl.harmonise_time(obs_list)

        if obsset == 'GHCNDEX_Rx1day':  # (varname_o=='Ann'):
            obs_raw = ghcnd_cat(obs_list)
        else:
            obs_raw = obs_list.concatenate_cube()  # Why do we have to take this step? We wouldn't if it was merge.
    else:
        obs_raw = iris.load_cube(obsname, nameconstraint & dateconstraint, callback=fl.datelabel_callback)
    fl.repair_coord_system(obs_raw)  # Repair issue with Iris 1.7 regridding
    if obsset.startswith('GPCC '):  # Fix lat units ("degrees_south") by replacing them with lon units
        obs_raw.coord('latitude').units = obs_raw.coord('longitude').units
        # We don't need to negate the coord points. I presume Iris has been clever without telling me.
        obs_raw.coord('longitude').circular = True  # This was set false, for some reason
    for axis in ['latitude', 'longitude']:
        if not obs_raw.coord(axis).has_bounds():
            obs_raw.coord(axis).guess_bounds()  # Do this before applying the region, else does not work for single grid boxes.
    obs = apply_region(region, obs_raw, varname_o)
    "Whether we used the constraint or the mask, these will propagate to the model through regridding and land masking."
    # Necessary for CRUTEM4, maybe not for other obs. Without these lines, area averaging won't work
    if (varname_o == 'tg') | (varname_o == 'pp'):  # EOBS temperatures are actuals
        return fl.anomalise(obs, clim_period[0], clim_period[1])
    else:
        return obs


def load_region_allforcings_ensemble(region, startyear, endyear, obs, modelroot, modelname, diagstr):
    '''Loads the regional all-forcings model data, concatenating RCP8.5 post-2005
    if this is necessary (the CMIP5 historical experiment finishes in 2005).'''
    freq = 'Amon'  # This should be passed in once we're working with dailies
    if (endyear > 2005) & (modelname != 'HadGEM3-A-N216'):
        if modelname in ['GISS-E2-H','GISS-E2-R']:#then we don't have enough RCP8.5. Use Ext instead
            print('Joining on historicalExt')
            histext_list = iris.cube.CubeList(
                [load_region_ensemble('historical', region, startyear, endyear, obs, modelroot, modelname, diagstr, freq),
                 load_region_ensemble('historicalExt', region, startyear, endyear, obs, modelroot, modelname, diagstr, freq)])
            all_ens = fl.harmonise_and_cat(histext_list)
        else:#Use RCP8.5
            print('Joining on RCP8.5')
            histext_list = iris.cube.CubeList(
                [load_region_ensemble('historical', region, startyear, endyear, obs, modelroot, modelname, diagstr, freq),
                 load_region_ensemble('rcp85', region, startyear, endyear, obs, modelroot, modelname, diagstr, freq)])
            all_ens = fl.harmonise_and_cat(histext_list)
    else:
        all_ens = load_region_ensemble('historical', region, startyear, endyear, obs, modelroot, modelname, diagstr,
                                       freq)
    return all_ens


def calc_bias(all_seas, obs_seas, startyear, eventyear, endyear):
    '''Calculates the mean offset between the model runs and the observations,
     over a climatology period running from the first year of data up to two years
     before the event we wish to attribute.'''
    if eventyear > endyear:
        lastclimyear = endyear
    else:
        lastclimyear = eventyear - 2
    obs_longclim = fl.calc_climatology(obs_seas, startyear, lastclimyear)
    # Produced the gridpoint-by-gridpoint full-period climatology for the obs. Now calc the mean difference from all-forcings.
    all_longclim = fl.calc_climatology(all_seas, startyear, lastclimyear).collapsed(['realization', 'physics_version'],
                                                                                    ia.MEAN)  # CUBE VERSION
    #    all_longclim=fl.calc_climatology(all_seas[0],startyear,eventyear-2)#CUBELIST
    #    for member in all_seas[1:]:
    #        all_longclim+=fl.calc_climatology(member,startyear,eventyear-2)
    # all_longclim/=len(all_seas)#Divide through to get ensemble mean
    all_longclim.units = 'kg m-2 month-1'  # IS THIS ALSO HAPPENING FOR tas?!
    for axis in ('time', 'month_number', 'year'):
        obs_longclim.coord(axis).attributes = all_longclim.coord(axis).attributes  # Fixes problem with GFDL model file
    bias = all_longclim.data - obs_longclim.data  # This won't work unless we regrid!
    return bias  # Removes time axis so we can subtract this from any timeseries


# def remove_bias(cube_ens,bias):
#    'Subtracts a given bias from each member of an ensemble (in a cube list).'
#    ens_bc=iris.cube.CubeList()
#    for member in cube_ens:
#        ens_bc.append(member-bias)
#    return ens_bc

# def event_area_mean(cube,eventyear):
#    yearconstraint=iris.Constraint(coord_values={'year': eventyear})
#    distn[pos:pos+n_membs[i]]=fl.area_mean(cube.extract(yearconstraint)).data.reshape(n_membs[i])

def event_ens_area_mean(ens_list, eventyear):
    '''Given a cubelist and a mask for an event, return an array of event means, 1 element per cube in the list.
    Note that this could be done as total for precip, but keeping as mean monthly over that season is just as good.'''
    # distn=np.empty(len(ens_list))#If there's likely to be members entirely missing data, use masked array
    n_membs = np.empty(len(ens_list))
    for i, cube in enumerate(ens_list):
        n_membs[i] = cube.coord('realization').points.size * cube.coord('physics_version').points.size
    distn = np.empty(n_membs.sum())
    pos = 0
    for i, cube in enumerate(ens_list):
        yearconstraint = iris.Constraint(coord_values={'year': eventyear})
        distn[pos:pos + n_membs[i]] = fl.area_mean(cube.extract(yearconstraint)).data.reshape(n_membs[i])
        pos += n_membs[i]
        # distn[i]=fl.area_mean(member[event]).data
    return distn


def get_ens_prob_wrt_thresh(distn,  # ens_list,#Now calculated outside this function
                            event, thresh, varn, operator='gt'):
    'Given a mask for the event year and a threshold, get the probability of exceeding the threshold in that year.'
    # distn=event_ens_area_mean(ens_list,event)
    if varn == 'pr':
        # Gamma version
        # shape,rate=fl.r_fit_gamma(distn)
        # prob=fl.r_gamma_lt_thresh(shape,rate,thresh[0])
        import scipy.stats as ss
        fit_alpha, fit_loc, fit_beta = ss.gamma.fit(distn, loc=0)
        prob = ss.gamma.cdf(thresh, fit_alpha, fit_loc, fit_beta)[0]
    elif varn == 'tas':
        # Gaussian part below. Modify this to include a switch from this to others, based on extremity?
        if type(thresh) == float:
            prob = fl.gaussian_lt_thresh(distn, thresh)
        else:
            prob = fl.gaussian_lt_thresh(distn, thresh[0])
        # prob=gev.pdf(thresh-mean, *gev.fit(distn))#Possible GEV example
    else:
        print("I don't know what distribution " + varn + " takes")

    #    if not 0.01<prob<0.99:#then this should be recalculated with generalised Pareto
    #        try:
    #            prob,_=fl.r_fit_covariate_gev(distn,thresh=0.99,event_value=thresh,type='GP')
    #        except ImportError:
    #            print("Couldn't find rpy2, using probability (%5.3f) from bulk distribution"%prob)

    if (operator == 'gt' or operator == '>' or operator == '>=' or operator == 'ge'):
        probout = 1 - prob  # 1-prob to give prob>obs (for heatwaves,floods etc)
    elif (operator == 'lt' or operator == '<' or operator == '<=' or operator == 'le'):
        probout = prob
    else:
        print(
            "I can only deal with probability greater than or less than the threshold. I assume you mean greater than.")
        probout = 1 - prob
    return probout  # ,distn


# def plot_forcing(anom,label,**kwargs):
#    from iris.plot import plot
#    plot(fl.area_mean(anom[0]).aggregated_by('year',ia.MEAN),label=label,**kwargs)
#    for i in range(1,anom.coord('physics_version').points.size):
#        plot(fl.area_mean(anom[i]).aggregated_by('year',ia.MEAN),**kwargs)


def plot_allvsnat_t_series(all_anom, nat_anom, obs, ylabel, lc=0,filename=None):
    'Plots each member of the All-forcings and Nat-forcings time series, plus the obs time series.'
    from matplotlib.pyplot import figure,close
    from pandas.plotting import register_matplotlib_converters
    register_matplotlib_converters()

    import iris.plot as ip
    fig = figure(num='NAT vs ALL ens members')
    subpl = fig.add_subplot(1, 1, 1)
    subpl.set_xlabel('Year')
    subpl.set_ylabel(ylabel)
    # For each of these loops, the first member is kept outside so that only one entry per forcing
    # is put in the legend. Nat-forcings plotted first so it's underneath the All-forcings.
    if lc == 0:
        labels = ['Without human influence', 'With human influence', 'Observed climate']  # labels=['Nat','All','Obs']
    else:
        labels = [None, None, None]
    #    if lc==0:
    # realizrange=range(all_anom.coord('realization').points.size)
    # physrange=range(all_anom.coord('physics_version').points.size)
    # print(realizrange,physrange)
    if type(all_anom) == iris.cube.CubeList:
        ip.plot(fl.area_mean(nat_anom[0][0]).aggregated_by('year', ia.MEAN), label=labels[0], color='#00b478')
        for nat_cube in nat_anom:
            for i in range(nat_cube.coord('physics_version').points.size):
                ip.plot(fl.area_mean(nat_cube[i]).aggregated_by('year', ia.MEAN), color='#00b478')
            # plot_forcing(nat_cube,labels[0],color='#00b478')
        ip.plot(fl.area_mean(all_anom[0][0]).aggregated_by('year', ia.MEAN), label=labels[1], color='#c80000')
        for all_cube in all_anom:
            # plot_forcing(all_cube,labels[1],color='#c80000')
            for i in range(all_cube.coord('physics_version').points.size):
                ip.plot(fl.area_mean(all_cube[i]).aggregated_by('year', ia.MEAN), color='#c80000')
    elif type(all_anom) == iris.cube.Cube:
        # plot_forcing(nat_anom,labels[0],color='#00b478')
        # plot_forcing(all_anom,labels[1],color='#c80000')
        ip.plot(fl.area_mean(nat_anom[0]).aggregated_by('year', ia.MEAN), label=labels[0], color='#00b478')

        # for i in range(nat_anom.coord('realization').points.size):
        for i in range(1, nat_anom.coord('physics_version').points.size):
            ip.plot(fl.area_mean(nat_anom[i]).aggregated_by('year', ia.MEAN), color='#00b478')
        #    if lc==0:
        ip.plot(fl.area_mean(all_anom[0]).aggregated_by('year', ia.MEAN), label=labels[1], color='#c80000')
        for i in range(1, all_anom.coord('physics_version').points.size):
            ip.plot(fl.area_mean(all_anom[i]).aggregated_by('year', ia.MEAN), color='#c80000')
    #    if lc==0:
    ip.plot(fl.area_mean(obs).aggregated_by('year', ia.MEAN), label=labels[2], color='#000000')
    #    else:
    #        ip.plot(fl.area_mean(obs).aggregated_by('year',ia.MEAN),color='#000000')
    subpl.legend(loc='upper left')
    if type(filename)==str:
        fig.savefig(filename)
        close(fig)
    else:
        ip.show(block=False)


def plot_allvsnat_distn(all_distn, obs_event, nat_distn=None, xlabel='T', varn='tas', comparator='gt',filename=None):
    '''Plots a Gaussian around each of the All and Nat distributions, marking each real point and the obs.'''
    n_sd_range = 4  # Number of standard deviations range of the distribution we'll plot.
    if type(nat_distn) == np.ma.core.MaskedArray:
        # Now take the max and min of both ranges together.
        Tranges = (np.array([[all_distn.mean()], [nat_distn.mean()]])
                   + np.array([-n_sd_range, n_sd_range]) * np.array([[all_distn.std()], [nat_distn.std()]]))
        # Written like this for simplicity. If we need to speed up that bit for large distns,
        # can pass mean and stdev from get_ens_prob_gt_thresh() instead of recalculating (and again in plot_distn).
        Trange = np.array([Tranges.min(), Tranges.max()])
    else:
        Trange = all_distn.mean() + np.array([-n_sd_range, n_sd_range]) * all_distn.std()
    # Generate an x axis for the plot from this
    Taxis = np.arange(Trange[0], Trange[1], (Trange[1] - Trange[0]) / 100.)
    import matplotlib.pyplot as mp
    fig = mp.figure(num='NAT vs ALL distributions')
    subpl = fig.add_subplot(1, 1, 1)
    subpl.set_ylabel('Probability density')
    subpl.set_xlabel(xlabel)
    allplot, allpeak = fl.plot_distn(Taxis, all_distn, varn, label='With human influence', obs=obs_event,
                                     color='#c80000',  # label='All'
                                     comparator=comparator, hatch='/')  # Plot All-forcings distribution
    if type(nat_distn) == np.ma.core.MaskedArray:
        natplot, natpeak = fl.plot_distn(Taxis, nat_distn, varn, label='Without human influence', obs=obs_event,
                                         color='#00b478',  # label='Nat'
                                         comparator=comparator, hatch='\\')  # then Nat-forcings
        dmax = max([allpeak, natpeak])
    else:
        dmax = allpeak
    mp.plot(np.array([obs_event, obs_event]), np.array([0, dmax]), label='Observed event',
            color='#000000')  # Vertical obs line#label='Obs'
    subpl.legend(loc='upper left')
    if type(filename)==str:
        fig.savefig(filename)
        mp.close(fig)
    else:
        mp.show(block=False)


def make_bootstrap_distn(distn, eventyear, obs_event, varn, comparator, n_bootstraps, sample_size):
    'Do a big bootstrap to get sample uncertainties on probabilities.'
    prob_distn = np.empty(n_bootstraps)
    for i in range(0, n_bootstraps):
        bs_distn = np.random.choice(distn, sample_size, replace=True)
        prob_distn[i] = get_ens_prob_wrt_thresh(bs_distn, eventyear,
                                                obs_event, varn, comparator)
    return prob_distn


def plot_bootstrap_distn(rr_distn, fractions=False,filename=None):
    'Takes a distribution of risk ratios and plots them as a log normal.'
    from matplotlib.pyplot import figure, show, close
    # from matplotlib.ticker import LogLocator,LogFormatter
    # from matplotlib.scale import Log2Transform,LogScale
    fig = figure(num='Risk ratio distribution')
    subpl = fig.add_subplot(1, 1, 1)
    subpl.set_ylabel('Probability density')
    subpl.set_xlabel('Risk Ratio')
    subpl.set_xscale('log', basex=2.0)
    # subpl.set_xscale(LogScale(Log2Transform()))
    # subpl.xaxis.set_transform(Log2Transform())
    # subpl.xaxis.set_major_locator(LogLocator(base=2.))
    # subpl.xaxis.set_major_formatter(LogFormatter(base=2.))
    n_sd_range = 4
    log_distn = np.log2(rr_distn)
    logrange = log_distn.mean() + np.array([-n_sd_range, n_sd_range]) * log_distn.std()  # What if mean or stdev are 0?
    # rrrange=rr_distn.mean()+np.array([-n_sd_range,n_sd_range])*rr_distn.std()
    # print(rrrange)
    # logrange=np.log2(rrrange)
    # rraxis=np.arange(rrrange[0],rrrange[1],(rrrange[1]-rrrange[0])/100.)
    print(logrange)
    print(2. ** logrange)
    rraxis = 2. ** (np.arange(logrange[0], logrange[1], (logrange[1] - logrange[0]) / 100.))
    rrplot, rrpeak = fl.plot_distn(rraxis, rr_distn, 'Risk Ratio', hatch='/',  # fill=False,
                                   label='Risk Ratio', color='#b8870a')  # Anthro brown
    if fractions:
        # xticks=subpl.xaxis.get_major_ticks()
        xticks = subpl.get_xticks()
        xlabels = np.vectorize(fl.reciprocal_label)(xticks)
        # subpl.xaxis.set_ticklabels(xlabels)
        subpl.set_xticklabels(xlabels)
    if type(filename)==str:
        fig.savefig(filename)
        close(fig)
    else:
        show(block=False)


def acplot(all_prob, nat_prob, obs_event, regionname, varn):
    '''Use Andy Ciavarella's routines to plot event bar charts.'''
    from event_def import Event as ACEvent  # Andy's event class
    from qfigs import barchart_p1p0  # Andy's simple barchart routine
    from matplotlib.pyplot import show
    ac_event = ACEvent([varn, 'season', regionname, obs_event],
                       nat_prob, all_prob)  # Include period logic later
    fig, ax = barchart_p1p0([ac_event], bounds=False)  # Include adjacent regions?
    show(block=False)


def loadandmask_cmip5da(startyear, endyear, region, obs, modelroot, modelname, varn, varname_o, freq):
    all_ens = load_region_allforcings_ensemble(region, startyear, endyear, obs[0], modelroot, modelname,
                                               varn)  # All forcings
    print('All-forcings ensemble loaded')
    nat_ens = load_region_ensemble('historicalNat', region, startyear, endyear, obs[0], modelroot, modelname, varn,
                                   freq)  # Natural forcings
    print('Natural-forcings ensemble loaded')
    # Below uses cubes preloaded into 15-member files. Faster on the desktop, slower on Pangeo
    # all_ens=load_and_regrid_cube('historical',region,startyear,endyear,obs[0],modelroot,modelname,varn,freq)#Doesn't include rcp85
    # nat_ens=load_and_regrid_cube('historicalNat',region,startyear,endyear,obs[0],modelroot,modelname,varn,freq)#Natural forcings

    notlandmask = fl.get_landsea_mask(modelroot + modelname + '/sftlf_fx_' + modelname + '_historical_r0i0p0.nc', obs)
    # notlandmask=get_landsea_mask(modelroot+'/sftlf_fx_'+modelname+'_historical_r0i0p0.nc',obs)
    mask_obs(obs, varn, varname_o, freq, notlandmask)
    # print('Masked observations')
    #all_clim = fl.calc_climatology(all_ens, clim_period[0], clim_period[1])#This uses each member's own climatology
    all_clim = fl.calc_climatology(all_ens, clim_period[0], clim_period[1]).collapsed(['realization','physics_version'],ia.MEAN)
    #The above (using the ens mean as climatology) should now work for different numbers of members to the climatology
    #from pdb import set_trace
    #set_trace()
    # print('Calculated climatology')
    all_anom = obs_and_landmask_anomalise_model(all_ens, varn, obs)#, clim=all_clim)
    # print('Masked all-forcings ensemble')
    nat_anom = obs_and_landmask_anomalise_model(nat_ens, varn, obs)#, clim=all_clim)
    # print('Masked natural-forcings ensemble')
    return all_anom, nat_anom  # ,Changes to obs masking should pass though automatically


def loadandmask_delayed(startyear, endyear, region, obs, modelroot, modelname, varn, varname_o, freq):
    from dask import delayed
    from dask.bag import from_delayed

    allload_delayed = delayed(load_region_allforcings_ensemble)(region, startyear, endyear, obs[0], modelroot,
                                                                modelname, varn)
    natload_delayed = delayed(load_region_ensemble)('historicalNat', region, startyear, endyear, obs[0], modelroot,
                                                    modelname, varn, freq)
    parallel_allnat_load = from_delayed([allload_delayed, natload_delayed])
    allnatlist = parallel_allnat_load.compute()
    all_ens, nat_ens = allnatlist

    notlandmask = fl.get_landsea_mask(modelroot + modelname + '/sftlf_fx_' + modelname + '_historical_r0i0p0.nc', obs)
    mask_obs(obs, varn, varname_o, freq, notlandmask)

    # all_clim=delayed(fl.calc_climatology)(allload_delayed,clim_period[0],clim_period[1])#Not previously done - make a keyword?
    all_clim = fl.calc_climatology(all_ens, clim_period[0], clim_period[1])  # Not previously done - make a keyword?
    # nat_anom=obs_and_landmask_anomalise_model(natload_delayed.compute(),varn,obs,clim=all_clim.compute())
    nat_anom = obs_and_landmask_anomalise_model(nat_ens, varn, obs, clim=all_clim)
    print('Natural-forcings ensemble loaded')
    all_anom = obs_and_landmask_anomalise_model(all_ens, varn, obs, clim=all_clim)
    # all_anom=obs_and_landmask_anomalise_model(allload_delayed.compute(),varn,obs,clim=all_clim.compute())
    print('All-forcings ensemble loaded')

    return all_anom, nat_anom


def yearpool(cube, startyear, endyear, years_to_pool=3):
    '''For using years before an event as additional ensemble members.
    Simplest method is to assume we've already done our annual (seasonal) means
    so we only have to copy elements either side to make multiple versions of the time axis
    then turn them into different realizations.'''
    from pdb import set_trace
    # tlen=len(cube.coord('time').points)
    # tdim=cube.coord_dims('time')[0]
    t_axis = cube.coord('time')[years_to_pool - 1:]
    year_axis = cube.coord('year')[years_to_pool - 1:]  # Could check this exists first
    enslen = len(cube.coord('realization').points)
    ens_axis = cube.coord('realization').copy()  # Definitely check for this if we plan to use this on obs
    ens_axis.points=np.arange(1,enslen+1)
    # We need to multiply realization up by physics version, then by time.
    # Which dimension is physics? Either use cube.coord_dims('physics_version') or use constraint in for loop.
    newcubelist = iris.cube.CubeList()
    for p in np.arange(0, cube.coord('physics_version').points.size) + 1:  #
        # Slice down physics dimension for 1 cube which is r,t,lat,lon
        # then append each to dimension. Do this for t as well. But is it a merge if len(r) is 1? Only first time?
        memberconstraint = iris.Constraint(coord_values={'physics_version': p})  # lambda phys: phys==p})
        i = 0
        while i < years_to_pool:
            # cubei=pslice[i:i+tlen-years_to_pool+1]#Offset data#THIS SHOULD BE DONE WITH CONSTRAINT TOO
            yearconstraint = iris.Constraint(
                coord_values={'year': lambda year: startyear + i <= year <= endyear - years_to_pool + i + 1})
            cubei = cube.extract(memberconstraint & yearconstraint)
            tdim = cubei.coord_dims('time')[0]  # Needs to be here, not before it's extracted, as we may lose dimensions in the extraction
            cubei.remove_coord('time')  # Put the years back to the original
            cubei.add_dim_coord(t_axis, tdim)  # Put the years back to the original#IS IT THIS AXIS?
            cubei.remove_coord('year')
            cubei.add_aux_coord(year_axis, tdim)
            if len(cubei.coords('realization')) > 0:
                cubei.remove_coord('realization')
            if len(cubei.coords('physics_version')) > 0:
                cubei.remove_coord('physics_version')

            newrealiz = ens_axis + enslen * i + enslen * years_to_pool * (p - 1)
            #This doesn't work without linear ensemble member numbers though, so the the change to ens_axis.points is done earlier

            if enslen > 1:  # there's already a realisation dimension
                cubei.add_dim_coord(newrealiz, 0)  # Increase ensemble member numbers for new "realization" number
            else:  # this was previously an ensemble of physics version, so we need to merge, not concatenate.
                cubei.add_aux_coord(newrealiz)
            # Consider using cube.coord_dims('realization' {or 'physics_version'}) if this hardwiring dimension is a problem
            newcubelist.append(cubei)  # Now append these "new" ensemble members
            i += 1
    #set_trace()
    if enslen > 1:
        return newcubelist.concatenate_cube()  # merge_cube()#Make into one cube with lots of members
    else:
        return newcubelist.merge_cube()


def yearpool_list(cubelist, years_to_pool=3):
    'Version of yearpool for ensembles constructed as lists instead of with a realization dimension.'
    tlen = len(cubelist[0].coord('time').points)
    t_axis = cubelist[0].coord('time')[years_to_pool - 1:]
    year_axis = cubelist[0].coord('year')[years_to_pool - 1:]  # Could check this exists first
    newcubelist = iris.cube.CubeList()
    i = 0
    while i < years_to_pool:
        for cube in cubelist:
            cubei = cube[i:i + tlen - years_to_pool + 1]  # Offset data
            cubei.remove_coord('time')  # Put the years back to the original
            cubei.add_dim_coord(t_axis, 0)  # Put the years back to the original
            cubei.remove_coord('year')
            cubei.add_aux_coord(year_axis, 0)
            newcubelist.append(cubei)
        i += 1
    return newcubelist


def var_inflation(correlation, obs, model_ens):  # ,model_mean,stdev_obs,stdev_ens,stdev_ensmean):
    'Prototype routine for variance inflation, based on Doblas-Reyes et al 2005 Tellus'
    stdev_obs = obs.aggregated_by('time', ia.STD_DEV)
    model_mean = np.mean(model_ens)  # Unless we have this elsewhere...
    stdev_model_mean = model_mean.aggregated_by('time', ia.STD_DEV)
    alpha = np.absolute(correlation) * stdev_obs / stdev_model_mean
    sbeta = np.sqrt(1 - correlation ^ 2) * stdev_obs
    model_infl = iris.cube.CubeList()
    for member in model_ens:
        stdev_member = member.aggregated_by('time', ia.STD_DEV)
        model_infl.append(alpha * model_mean + sbeta / stdev_member * (member - model_mean))
    return model_infl


def smoothed_globalobs(startyear, endyear, varo, obsset, startmonth, endmonth, aggregator):
    'This is for using as a covariate with van Oldenborgh or Jeon-style extreme value fits.'
    globalobs = load_obs(startyear, endyear, False, obsset, varo)
    # globalobs_seas=fl.filter_months(globalobs,startmonth,endmonth).aggregated_by('year',aggregator)#Should we subselect seasonal?
    globalobs_ann = globalobs.aggregated_by('year', aggregator)  # Should we subselect seasonal?
    globalavg_obs = fl.area_mean(globalobs_ann)
    # GJvO does a 3-year running mean on T_global. We will too (but Jeon does the 13-point filter from Solomon et al 2007).
    return globalavg_obs.rolling_window('time', aggregator, 3)


def load_process_data(startyear, endyear, eventyear, startmonth, endmonth, region, obs, modelroot, modelname, varn,
                      varname_o, freq, lc):
    #    if iris.__version__.startswith('2'):
    ##        from dask import delayed
    ##        loadandmask_delayed=delayed(loadandmask_cmip5da)
    #        anomalies=loadandmask_delayed(startyear,endyear,region,obs,modelroot,modelname,varn,varname_o,freq)
    #        all_anom,nat_anom=anomalies#.compute()#This should now run loadandmask_cmip5da in parallel
    #    else:
    all_anom, nat_anom = loadandmask_cmip5da(startyear, endyear, region, obs, modelroot, modelname, varn, varname_o,
                                             freq)  # GENERIC ATTRIBUTION load routine

    # Look at a single season. Event attribution example. Make sure this is within the model run period!
    obs_seas = fl.filter_months(obs, startmonth, endmonth).aggregated_by('year',ia.MEAN)  # Mean for each year over the season of interest
    all_seas = fl.seasonal_ensemble(all_anom, startmonth, endmonth, obs_seas.coord('time'))
    nat_seas = fl.seasonal_ensemble(nat_anom, startmonth, endmonth, obs_seas.coord('time'))

    # Apply bias correction from the years before the one we're interested in
    #bias = calc_bias(all_seas, obs_seas, startyear, eventyear, endyear)
    bias = calc_bias(all_seas, obs_seas, clim_period[0], eventyear, clim_period[1])
    all_seas_bc = all_seas - bias
    nat_seas_bc = nat_seas - bias
    print('Using climatology period to correct bias')
    #print('Ignoring the bias')

    # Variance inflation could go here, once we've developed it for attribution
    # Calculate obs correlation with All for pre- and post-change. Also obs std dev and All pre/post?
    # (Omar: ensemble mean Pearson correlation in time (interannual))
    # Are any of the variances from detrended data?
    # Apply pre- to Nat ensemble, post- to All.

    # Produce distributions for the event year
    event = (obs_seas.coord('year').points == eventyear)
    obs_event = fl.area_mean(obs_seas[event]).data
    return all_seas_bc, nat_seas_bc, event, obs_event, bias, all_seas, nat_seas, obs_seas


def load_process_eucleia_short(startyear, endyear, eventyear, startmonth, endmonth, region, obs, modelroot, modelname,
                               varn, varname_o, freq, lc):
    notlandmask = fl.get_landsea_mask(modelroot + modelname + '/sftlf_fx_' + modelname + '_historical_r0i0p0.nc', obs)
    mask_obs(obs, varn, varname_o, freq, notlandmask)
    # This is a tricky one. We now need obs to encompass the eventyear.
    # If that's beyond the endyear of the climatology, we need a truncated version of the obs here.
    obs_long = obs[obs.coord('year').points <= endyear]
    obs_short = obs[obs.coord('year').points == eventyear]
    all_long = load_region_allforcings_ensemble(region, startyear, endyear, obs_long[0], modelroot, modelname,
                                                varn)  # All forcings
    all_long_anom = obs_and_landmask_anomalise_model(all_long, varn, obs_long)  # Need to get climatology from here
    all_clim = fl.calc_climatology(all_long_anom, clim_period[0], clim_period[1]).collapsed(['physics_version'],
                                                                                            ia.MEAN)
    if (eventyear == 2014 or eventyear == 2015):
        all_label = 'historicalShort'
        nat_label = 'historicalNatShort'
    elif eventyear == 2016 and endmonth <= 5:  # Ext01 was done for JFMAM2016
        all_label = 'historicalExt01'
        nat_label = 'historicalNatExt01'
    else:
        all_label = 'historicalExt'
        nat_label = 'historicalNatExt'
    all_short = load_region_ensemble(all_label, region, eventyear, eventyear, obs[0], modelroot, modelname, varn,
                                     freq, startendmonths=[startmonth,endmonth])  # All forcings - needs clim
    #    print(obs_short)
    #    print(obs_short.coord('month_number'))
    #    print(obs_short.coord('day_of_year'))
    #    print(all_short)
    #    print(all_short.coord('month_number'))
    #    print(all_short.coord('day_of_year'))
    #    print(all_short.coord('month_number').points.min())
    #    print(all_short.coord('month_number').points.max())
    obs_short_cropped = fl.filter_months(obs_short, int(all_short.coord('month_number').points.min()),
                                         int(all_short.coord('month_number').points.max()), yearshift=False)
    all_clim_cropped = fl.filter_months(all_clim, int(all_short.coord('month_number').points.min()),
                                        int(all_short.coord('month_number').points.max()), yearshift=False)
    all_short_anom = obs_and_landmask_anomalise_model(all_short, varn, obs_short_cropped, clim=all_clim_cropped)
    nat_short = load_region_ensemble(nat_label, region, eventyear, eventyear, obs[0], modelroot, modelname, varn, freq)
    nat_short_anom = obs_and_landmask_anomalise_model(nat_short, varn, obs_short_cropped, clim=all_clim_cropped)

    obs_seas = fl.filter_months(obs, startmonth, endmonth).aggregated_by('year',
                                                                         ia.MEAN)  # Mean for each year over the season of interest
    obs_long_seas = obs_seas[obs_seas.coord('year').points <= endyear]
    obs_short_seas = obs_seas[obs_seas.coord('year').points == eventyear]
    all_long_seas = fl.seasonal_ensemble(all_long_anom, startmonth, endmonth, obs_long_seas.coord('time'))
    all_short_seas = fl.seasonal_ensemble(all_short_anom, startmonth, endmonth, obs_short_seas.coord('time'))
    nat_short_seas = fl.seasonal_ensemble(nat_short_anom, startmonth, endmonth, obs_short_seas.coord('time'))

    bias = calc_bias(all_long_seas, obs_seas, startyear, eventyear,
                     endyear)  # This is a spatially varying bias. Is that what we want?
    all_seas_bc = all_short_seas - bias
    nat_seas_bc = nat_short_seas - bias

    # xlabel=all_long_anom[0].long_name
    # plot_allvsnat_t_series(all_seas_bc,nat_seas_bc,obs_seas,xlabel,lc)#Time series plots
    event = (obs_seas.coord('year').points == eventyear)
    obs_event = fl.area_mean(obs_seas[event]).data
    return all_seas_bc, nat_seas_bc, event, obs_event, bias, all_short_seas, nat_short_seas, obs_seas


# =====================
# MAIN CODE STARTS HERE
# =====================
def main(modelname):
    from time import time
    calcstart = time()

    obs_event = None  # Read it from obs time series, unless pre-defined below

    #    from matplotlib import use
    #    use('Agg')#Stop Matplotlib from using Xwindows

    #    startyear=1960
    #    endyear=2003#10
    #    #The following timings are now around 7x faster due to optimising with callbacks and constraints within load routines
    #    region=fl.giorgi_coord(15)#with CSIRO, took 4 mins 30 for region 15, 3:20 for 11
    #    #region=fl.get_srex_mask(13)#Runs out of memory with this method
    #    #region=fl.get_srex_shape(17)#Took around 3 mins 20 for region 13
    #    #print 'Examining '+region.attributes['name']+' (SREX region %2d)'%region.attributes['number']
    #    eventyear=2002
    #    startmonth=6
    #    endmonth=8
    #    #eventyear=1999
    #    #startmonth=12
    #    #endmonth=2
    #    varn='pr' # Precipitation='pr', Temperature='tas'
    #    pool=True

    # Summer 2012 N Europe (e.g. UK) flood
    #    startyear=1960
    #    endyear=2012
    #    eventyear=2012
    #    startmonth=6
    #    endmonth=8
    #    varn='pr'
    #    region=fl.get_srex_shape(11)
    #    print 'Examining '+region.attributes['name']+' (SREX region %2d)'%region.attributes['number']
    #    pool=True

    #    #Setup for Summer 2012 European heatwave
    #    startyear=1960
    #    endyear=2012
    #    startmonth=6
    #    endmonth=8
    #    thresh=0.333#for pp and rr#0.667#for tg#
    #    varo='rr'#'tg'#'pp'#
    #    obsset='E-OBS'
    #    #region=[-10,50,25,60]#northern Europe
    #    region=[-10,35,25,45]#southern Europe = -10 to 25E, 35 to 45N
    #    #region=[-10,35,25,60]#all Europe
    #    eventyear=2012
    #    pool=True

    #    #Summer 2015 C Europe heatwave
    #    startyear=1960#for tas#1984#for mrso#
    #    endyear=2013
    #    startmonth=6
    #    endmonth=8
    #    varo='tg'#'pp'#'rr'#
    #    thresh=0.667#for tg and pp#0.333#for rr and SWBM#
    #    obsset='E-OBS'#'SWBM'#
    #    #region=False#Global average
    #    region=fl.get_srex_shape(12).union(fl.get_srex_shape(13))
    #    print 'Examining '+region.attributes['name']+' (SREX region '+region.attributes['number']+')'
    #    eventyear=2015
    #    pool=False

    # Warm extremes for defrosting reindeer. Summer 2016, Siberia SREX region but north of Noyabr'sk
    #    startyear=1960
    #    endyear=2013
    #    startmonth=6
    #    endmonth=8
    #    thresh=0.667
    #    obsset='CRUTEM4'
    #    varo='temperature_anomaly'
    #    region=[40,66,180,77.5]#North of Salekhard (closest to outbreak) and south of Cape Chelyuskin (northernmost mainland point).
    #    eventyear=2011#2016#
    #    pool=False
    # Further tests:
    # region=[-10,50,25,60]#northern Europe #Check region that crosses GMT works
    # thresh=0.333
    # obsset='GPCC 2.5 degree'
    # varo='precip'

    # China 2017 heatwave?
    # startyear=1960
    # endyear=2012
    # region=[106,25,122,37]#[25-37N, 106-122E]
    # Pentad Tmax - not yet implemented

    # UK 2003 or 2010 heatwaves, SSTs most analoguous to 2018 and 2019
    startyear = 1960
    endyear = 2012
    #startmonth = 3
    #endmonth = 5
    startmonth=6
    endmonth=6
    thresh = 0.667
    obsset = 'CRUTEM4'
    varo = 'temperature_anomaly'
    region=fl.giorgi_coord(11)#Mediterranean
    regionname='Mediterranean Basin'
    #region = fl.giorgi_coord(12)
    #regionname = 'Northern Europe'
    #    obsset='E-OBS'
    #    varo='tg'
    # obs_event=3.05#2018 UK provisional anom from NCIC is 3.2 in May and 2.9 in June
    # obsset='GPCC 1.0 degree'
    # varo='precip'
    # thresh=0.333
    # region=[-6,50,2,59]#GB, excludes Hebrides and Ireland. Doesn't overlap well with CRUTEM grid though.
    # region=[-4.99,50.01,-0.01,59.99]#Prevents interaction with bounds so we just get 2 CRUTEM grid boxes.
    eventyear =2019#2003#2018# 2010  #
    pool = False#True  #

    # Now translate the obsset into filename and CF variable names. Should this next bit be done with an object instead?
    obsfname, varo, varn, freq = fl.set_variables(obsset, varo)

    if endyear > eventyear:
        obs = load_obs(startyear, endyear, region, obsset, varo)  # ,obsfname)#Observations (already anomalies)
    else:
        obs = load_obs(startyear, eventyear, region, obsset, varo)
    multi_all_seas_bc = iris.cube.CubeList()  # Multi-model ensembles
    multi_nat_seas_bc = iris.cube.CubeList()

    for i, modeln in enumerate(modelname):
        print('Loading model ' + modeln)
        if modeln == 'HadGEM3-A-N216':
            # modelroot='/data/local/hadlf/'+modeln+'/'
            # modelroot='/project/detn/'+modeln+'/historical/'+varn+'/'+freq+'/'#Except we need "historical" to change too...
            modelroot = '/project/detn/'
            # modelroot='/data/cr1/hadlf/'
        else:
            #modelroot = '/data/local/hadlf/data_els028/CMIP5/'  # +modeln+'/'
            modelroot = '/data/users/hadlf/CMIP5/'  # +modeln+'/'
        if endyear >= eventyear:
            modeldata = load_process_data(startyear, endyear, eventyear, startmonth, endmonth, region, obs, modelroot,
                                          modeln, varn, varo, freq, i)
        else:
            modeldata = load_process_eucleia_short(startyear, endyear, eventyear, startmonth, endmonth, region, obs,
                                                   modelroot, modeln, varn, varo, freq, i)
        all_seas_bc, nat_seas_bc, event, obs_event_read, bias, all_seas, nat_seas, obs_seas = modeldata
        multi_all_seas_bc.append(all_seas_bc)  # was extend for list
        multi_nat_seas_bc.append(nat_seas_bc)

        if endyear >= eventyear:
            xlabel = varo  # all_seas_bc[0].long_name
            plot_allvsnat_t_series(all_seas_bc, nat_seas_bc, obs_seas, xlabel,
                                   i)  # Time series plots - move outside loop?

    if obs_event == None:  # then it hasn't been manually defined.
        obs_event = obs_event_read  # Use the obs value (and threshold) from the obs time series

    duration = time() - calcstart
    print('Read time %3d minutes and %2d seconds' % (duration // 60, duration % 60))

    'Year pooling a la King et al. Needed some work to make it compatible with Jeon et al quantile-matching.'
    if pool:
        years_to_pool = 10
        all_pooled = iris.cube.CubeList()
        for model in multi_all_seas_bc:
            all_pooled.append(yearpool(model, startyear, endyear, years_to_pool))
        nat_pooled = iris.cube.CubeList()
        for model in multi_nat_seas_bc:
            nat_pooled.append(yearpool(model, startyear, endyear, years_to_pool))
    else:
        years_to_pool = 1
        all_pooled = multi_all_seas_bc
        nat_pooled = multi_nat_seas_bc

    #for cube in all_pooled:
    #    print(cube.coord('realization'))
    #for coord in ['latitude','longitude','realization','time']:
    #    print(coord)
    #    fl.coord_diffs(all_pooled[0].coord(coord), all_pooled[1].coord(coord))
    #from iris.util import describe_diff
    #describe_diff(all_pooled[0],all_pooled[1])
    #print(fl.cubise_ensemble(all_pooled).coord('realization').points)
    all_distn = fl.area_mean(fl.cubise_ensemble(all_pooled)).data[:, :-1]  # Don't include final year
    # distn_shape=all_distn.shape
    all_distn = all_distn.reshape(all_distn.size)  # Flatten the distribution for R
    nat_distn = fl.area_mean(fl.cubise_ensemble(nat_pooled)).data[:, :-1]  # Don't include final year
    nat_distn = nat_distn.reshape(nat_distn.size)  # Flatten

    # if pool:
    #     try:
    #         "Jeon et al (with or without the pooling?)."
    #         globalrunning3 = smoothed_globalobs(startyear, endyear, varo, obsset, startmonth, endmonth, ia.MEAN)
    #         global_on_eventyear = globalrunning3[globalrunning3.coord('year').points == eventyear - 2]
    #         obsprob_cov, obsfit = fl.r_fit_covariate_gev(fl.area_mean(obs_seas).data[1:-1], thresh, globalrunning3.data,
    #                                                      event_value=obs_event,
    #                                                      event_covariate=global_on_eventyear.data)
    #         # We have to drop the first two years to match the data to the 3-year running average
    #         # Bring ensemble members into a cube dimension so they can easily form part of the distrubution
    #         globalrepeat = np.tile(globalrunning3.data[years_to_pool - 2:], all_distn.size / (
    #                     endyear - startyear - years_to_pool + 1))  # Repeat the global series for covariate location
    #         all_returnlevel, allfit = fl.r_fit_covariate_gev(all_distn, thresh, globalrepeat, obs_prob=obsprob_cov,
    #                                                          event_covariate=global_on_eventyear.data)  # Return level for the observed probability
    #         if all_returnlevel != None:
    #             natprob_cov, natfit = fl.r_fit_covariate_gev(nat_distn, thresh,
    #                                                          event_value=all_returnlevel)  # Natural prob for the all-forcings level
    #             print('Jeon et al method P_obs = %5.3f' % obsprob_cov, 'P_nat = %5.3f' % natprob_cov)
    #             dblp_cov = np.log2(obsprob_cov / natprob_cov)
    #             print('DBLP by Jeon method = %5.2f' % dblp_cov)
    #             print()
    #         else:
    #             natfit = None
    #     except ModuleNotFoundError:
    #         print("rpy2 is missing, so couldn't use R for Jeon method")
    obsfit = allfit = natfit = None
    global_on_eventyear = None


    # When we have all that working, we need to work out the uncertainty too!

    # Save all data
    #    iris.save(all_ant_seas_bc, varn+'_'+str(region)+'_'+str(eventyear)+'_all_ant_members.nc')
    #    iris.save(all_nat_seas_bc, varn+'_'+str(region)+'_'+str(eventyear)+'_all_nat_members.nc')
    #    iris.save(obs_seas,varn+'_'+str(region)+'_'+str(eventyear)+'_obs.nc')

    'Thresh, defined above in the code, is the outer threshold for GEV etc, and for now'
    "the way of deciding whether we're interested in greater or less than the obs event."
    if thresh < 0.5:
        comparator = 'lt'
    else:
        comparator = 'gt'
    # We now calculate *_distn above
    all_prob = get_ens_prob_wrt_thresh(all_distn, eventyear, obs_event, varn, comparator)  # event[years_to_pool-1:]
    nat_prob = get_ens_prob_wrt_thresh(nat_distn, eventyear, obs_event, varn, comparator)  # event[years_to_pool-1:]

    xlabel = varo  # all_seas_bc[0].long_name
    plot_allvsnat_distn(all_distn,obs_event,nat_distn,xlabel,varn,comparator)

    #plot_allvsnat_distn(all_distn, obs_event, None, varo, varn, comparator)  # Just all-forcings

    if all_prob == 0:
        if nat_prob == 0:
            rr = dblp = far = 0
        else:
            far = -np.inf
            rr = 0
            dblp = -np.inf
    else:
        if nat_prob == 0:
            far = 1
            rr = np.inf
            dblp = np.inf
        else:
            far = 1 - (nat_prob / all_prob)  # Fraction of attributable risk
            rr = all_prob / nat_prob
            dblp = np.log2(rr)
    psufficient = 1 - (1 - all_prob) / (1 - nat_prob)
    print('P_all = %5.3f' % all_prob, 'P_nat = %5.3f' % nat_prob)
    print('Fraction of Attributable Risk (FAR) = %5.3f' % far)
    print('Risk ratio = %5.1f' % rr)
    print('Difference of binary logarithms of probability (DBLP) = %5.2f' % dblp)
    print('Probability of sufficient causality = %5.3f' % psufficient)
    duration = time() - calcstart
    print('Calculation time %3d minutes and %2d seconds' % (duration // 60, duration % 60))

    # Do a big bootstrap to get sample uncertainties on probabilities
    n_bootstraps = 500
    sample_size = 525
    all_prob_distn = make_bootstrap_distn(all_distn, eventyear, obs_event, varn, comparator, n_bootstraps, sample_size)
    nat_prob_distn = make_bootstrap_distn(nat_distn, eventyear, obs_event, varn, comparator, n_bootstraps, sample_size)
    rr_distn = all_prob_distn / nat_prob_distn
    plot_bootstrap_distn(rr_distn, fractions=True)

    acplot(all_prob, nat_prob, obs_event, regionname, varn)

    if pool:
        return obsfit, allfit, natfit, global_on_eventyear
    else:
        return None, None, None, None


#    np.savetxt(varn+'_'+str(region)+'_'+str(eventyear)+'_all_ant_members.txt', np.column_stack((all_distn)))
#    np.savetxt(varn+'_'+str(region)+'_'+str(eventyear)+'_all_nat_members.txt', np.column_stack((nat_distn)))
#    np.savetxt(varn+'_'+str(region)+'_'+str(eventyear)+'_obs.txt', np.column_stack((obs_event)))

# WHEN CODE FIRST RUNS, IT DOES THE IMPORTS THEN SKIPS PAST ALL THE FUNCTIONS DOWN TO HERE
# Code below only determines models to use based on being run from commandline with a list of model names.
# It then calls the main code above.
if __name__ == '__main__':
    from sys import argv

    modelnames = [''] * (len(argv) - 1)
    if len(argv) > 1:
        for a in range(1, len(argv)):
            modelnames[a - 1] = argv[a]
    else:
        modelnames = ['HadGEM3-A-N216']  # ['HadGEM2-ES','CSIRO-Mk3-6-0','CanESM2']##
    obsfit, allfit, natfit, global_on_eventyear = main(modelnames)
# from rpy2.robjects.packages import importr
# extRemes=importr('extRemes')
# allparams=extRemes.findpars(allfit)
# print allparams
