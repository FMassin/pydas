import numpy as np
import h5py
from obspy import UTCDateTime, Stream, Trace
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.clients.nrl import NRL
import string
from  pandas import ExcelFile
from os import system
from glob import glob

# default seiscomp schema for h5toseed
defaultseiscompschema = 'xmlns="http://geofon.gfz-potsdam.de/ns/seiscomp3-schema/0.12" version="0.12"'

# default binding template for h5toseed
defaultbindingtemplate = '''      <module publicID="Config/trunk" name="trunk" enabled="true">
                <station publicID="Config/trunk/%s/%s"
                            networkCode="%s" stationCode="%s" enabled="true">
                    <setup name="default" enabled="true">
                    <parameterSetID>ParameterSet/trunk/Station/%s/%s/default</parameterSetID>
                    </setup>
                </station>
            </module>
            <parameterSet publicID="ParameterSet/trunk/Station/%s/%s/default" created="1970-01-01T00:00:00.000000Z">
                <moduleID>Config/trunk</moduleID>
                <parameter publicID="smi:ch.ethz.sed/h5toseed/Parameter/19700101000000.000000.00000">
                    <name>detecLocid</name>
                    <value>"00"</value>
                </parameter>
                <parameter publicID="smi:ch.ethz.sed/h5toseed/Parameter/19700101000000.000000.00001">
                    <name>detecStream</name>
                    <value>%s%s%s</value>
                </parameter>
            </parameterSet>
        '''

def scexec(scmodule='scolv',
           d='postgresql://***:***@eq20a.ethz.ch:5432/sc3dba?column_prefix=m_',
           configdb='conf.merged.xml',
           inventorydb='inv.merged.xml',
           recordstream='routing://file/data.mseed??match=OF.*.*.*;fdsnws/arclink.ethz.ch:8080/fdsnws/dataselect/1/query??match=*.*.*.*',
           options="-u test -E 'smi:ch.ethz.sed/sc3a/2019qofceu' --debug",
           ssh=None):
    """
    Executes the SeisComP module command with the given parameters.

    Args:
        scmodule (str, optional): The name of the SeisComP module to execute. Defaults to 'scolv'.
        d (str, optional): The database connection string. Defaults to 'postgresql://***:***@eq20a.ethz.ch:5432/sc3dba?column_prefix=m_'.
        configdb (str, optional): The path to the merged configuration XML file. Defaults to 'conf.merged.xml'.
        inventorydb (str, optional): The path to the merged inventory XML file. Defaults to 'inv.merged.xml'.
        recordstream (str, optional): The recordstream URL for accessing data. In order to get OF data from file and the rest from fdsnws, the defaults is set to 'routing://file/data.mseed??match=OF.*.*.*;fdsnws/arclink.ethz.ch:8080/fdsnws/dataselect/1/query??match=*.*.*.*'.
        options (str, optional): Additional options for the SeisComP module command. Defaults to "-u test -E 'smi:ch.ethz.sed/sc3a/2019qofceu' --debug".

    Returns:
        None
    """
    if ssh is None:
        system(f"{scmodule} --offline -d '{d}'  --config-db file://{configdb}  --inventory-db file://{inventorydb} --recordstream='{recordstream}' {options}")
    ###MAKE AN SSH VERSION?
    ##!ssh sc20ag.ethz.ch scolv -u test  -d 'postgresql://***:***@eq20a.ethz.ch:5432/sc3dba?column_prefix=m_'  --config-db file:///home/fmassin/optic-fiber/data/FLUELA/conf.merged.xml  --inventory-db file:///home/fmassin/optic-fiber/data/FLUELA/inv.merged.xml --recordstream='routing://file//home/fmassin/optic-fiber/data/FLUELA/data.mseed??match=OF.*.*.*\;sdsarchive//rz_nas/miniseed/??match=*.*.*.*' -E 'smi:ch.ethz.sed/sc20a/Event/2022dxmayg' --debug

def seedtosc(metadatastationxml='metadata.stationxml',
             metadatascxml='metadata.scxml',
             db2mergewith='postgresql://****:****@eq20a.ethz.ch:5432/sc3dba?column_prefix=m_',
             confxml='conf.xml',
             invmergedxml='inv.merged.xml',
             confmergedxml='conf.merged.xml'):
    """
    Convert metadata from FDSNXML format to SCXML format and merge it within SCXML files.

    Notes:
        SeisComP utility `scxmldump` needs to be available and properly configured on your system

    Args:
        metadatastationxml (str): Path to the FDSNXML metadata file.
        metadatascxml (str): Path to the output SCXML metadata file.
        db2mergewith (str): Database connection string for merging with existing SCXML files.
        confxml (str): Path to the configuration XML file.
        invmergedxml (str): Path to the merged inventory XML file.
        confmergedxml (str): Path to the merged configuration XML file.
    """

    system(f"fdsnxml2inv '{metadatastationxml}'  '{metadatascxml}'")
    system(f"scxmldump -d '{db2mergewith}' -I | scxmlmerge - '{metadatascxml}' > '{invmergedxml}'")

    if confxml is None:
        return
    system(f"scxmldump -d '{db2mergewith}' -C | scxmlmerge - '{confxml}'       > '{confmergedxml}'")

    ### MAKE AN SSH VERSION?
    ##!rsync -avzl ../../../optic-fiber/ sc20ag.ethz.ch:optic-fiber/ 
    ##!ssh sc20ag.ethz.ch fdsnxml2inv optic-fiber/data/FLUELA/metadata.stationxml  optic-fiber/data/FLUELA/metadata.scxml
    ##!ssh sc20ag.ethz.ch scxmldump -d 'postgresql://****:****@eq20a.ethz.ch:5432/sc3dba?column_prefix=m_'  -I \| scxmlmerge - optic-fiber/data/FLUELA/metadata.scxml  \> optic-fiber/data/FLUELA/inv.merged.xml
    ##!ssh sc20ag.ethz.ch scxmldump -d 'postgresql://****:****@eq20a.ethz.ch:5432/sc3dba?column_prefix=m_'  -C \| scxmlmerge - optic-fiber/data/FLUELA/conf.xml        \> optic-fiber/data/FLUELA/conf.merged.xml
    


def h5toseed(strainfile=None,
             velocityfile=None,
             coordinates=None,
             locationsperstation = 4,
             sensitivity=1E7,
             orientation_code='F',
             strain_code='S', # linear strainmeter  https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/#linear-strain and https://www.unavco.org/data/strain-seismic/bsm-data/lib/docs/BOREHOLE_StrainSheet_DP_11.pdf
             band_code='H',
             confxml="conf.xml",
             datamseed='data.mseed',
             metadatastationxml="metadata.stationxml",
             nrlbackup=None,
             datakey='/data',
             stationdim=1,
             bindingtemplate=defaultbindingtemplate,
             seiscompschema=defaultseiscompschema,
             **attr):
    """
    Convert HDF5 files to SEED format.

    Args:
        strainfile (str): Path to the HDF5 file containing strain data.
        velocityfile (str): Path to the HDF5 file containing velocity data.
        coordinatefile (str): Path to the Excel file containing coordinate information.
        locationsperstation (int): Maximum number of locations per station.
        sensitivity (float): Sensitivity value for scaling the data.
        orientation_code (str): Orientation code for the SEED channel naming.
        strain_code (str): Strain code for the SEED channel naming.
        band_code (str): Band code for the SEED channel naming.
        confxml (str): Path to the output SeisComP station configuration XML file. Skip SeisComP bindings by setting to None.
        datamseed (str): Path to the output MiniSEED file.
        metadatastationxml (str): Path to the output StationXML file.
        nrlbackup (str): Path to the NRL backup file.
        datakey (str): Key to access the data in the HDF5 file.
        bindingtemplate (str): Template for generating SeisComP station configuration bindings. 
        seiscompschema (str): SeisComP schema version. 
        **attr: Additional attributes for the SEED metadata.

    Returns:
        tuple: A tuple containing the Inventory object, Stream object, and the SeisComP station configuration XML string.
    """

    alphanum = [l+k for l in '0123456789' for k in '0123456789']
    alphanum += [l+k for l in string.ascii_uppercase for k in string.ascii_uppercase] 
    alphanum += [l+k for l in string.ascii_lowercase for k in string.ascii_lowercase]     
    print('max num of stations:',len(alphanum))

    # Create all the base objects. 
    attr['starttime'] = UTCDateTime(attr['starttime'])

    if coordinates is not None:
        try:
            xls = ExcelFile(r"%s"%coordinates)
            sheetX = xls.parse(0)
            latitudes = sheetX['Var1']
            longitudes = sheetX['Var2']
            elevations = sheetX['Var3']
        except:
            pass
        try:
            latitudes = coordinates[0]
            longitudes = coordinates[1]
        except:
            pass
        try:
            elevations = coordinates[2]
        except:
            elevations = [attr['depth'] for l in latitudes]

    stream = Stream()

    # These strongly follow the hierarchy of StationXML files.
    inv = Inventory(
        # We'll add networks later.
        networks=[],
        # The source should be the id whoever create the file.
        source=attr['source'])

    net = Network(
        # This is the network code according to the SEED standard.
        code=attr['network'],
        # A list of stations. We'll add one later.
        stations=[],
        description=attr['description'],
        # Start and end dates are optional.
        start_date=attr['starttime'])

    sta = Station(
        # This is the station code according to the SEED standard.
        code=attr['station'],
        latitude=attr['latitude'],
        longitude=attr['longitude'],
        elevation=attr['elevation'],
        creation_date=UTCDateTime(2016, 1, 2),
        site=Site(name=attr['description']))


    # By default this accesses the NRL online. Offline copies of the NRL can
    # also be used instead
    # The contents of the NRL can be explored interactively in a Python prompt,
    # see API documentation of NRL submodule:
    # http://docs.obspy.org/packages/obspy.clients.nrl.html
    # Here we assume that the end point of data logger and sensor are already
    # known:
    velocity_response = NRL(nrlbackup).get_response(sensor_keys=['Generic', 'Unity Velocity Sensor'],
                                datalogger_keys=['Generic', 'Unity']) 
    #https://ds.iris.edu/NRL/sensors/generic/RESP.XX.NS000..BHZ.UNITY.DC.1

    strain_response = NRL(nrlbackup).get_response(sensor_keys=['Generic', 'Unity Velocity Sensor'],
                                datalogger_keys=['Generic', 'Unity']) 
    #https://ds.iris.edu/NRL/sensors/generic/RESP.XX.NS000..BHZ.UNITY.DC.1

    strain_response.instrument_sensitivity.input_units = attr['units'] 
    strain_response.instrument_sensitivity.input_units_description = attr['units_description'] 
    strain_response.response_stages[0].input_units = attr['units'] 
    strain_response.response_stages[0].input_units_description = attr['units_description']


    # Read the HDF5 files
    if strainfile is not None:
        strain_array = None
        for f in np.sort(glob(strainfile)):
            print(f)
            with h5py.File(f, 'r') as file:
                tmpdata = np.array(file[datakey])
                if stationdim ==1:
                    tmpdata = np.transpose(np.array(file[datakey]))

                if strain_array is None:
                       strain_array = tmpdata / sensitivity
                else:
                    strain_array = np.concatenate([strain_array, tmpdata / sensitivity], axis=0)

                print('strain_array.shape : %s'%str(strain_array.shape))
        
        strain_array = strain_array.astype(np.float64)

    if velocityfile is not None:
        with h5py.File(velocityfile, 'r') as file:
            velocity_array = np.transpose(np.array(file[datakey])) / sensitivity
            print('velocity_array.shape : %s'%str(velocity_array.shape))

        velocity_array = velocity_array.astype(np.float64)

    stations = []
    locationindex = -1
    stationindex = -1
    sta = None
    
    for i in range(strain_array.shape[1]): 
        
        locationindex += 1

        # Create new station once max number of location codes reached
        if locationindex >= locationsperstation or stationindex == -1:
            
            if sta is not None:
                net.stations.append(sta)

            locationindex = 0
            stationindex += 1

            attr['station'] = attr['station'][:3]+alphanum[stationindex]
            stations += [attr['station']]
        
            # This is the station code according to the SEED standard.
            sta = Station(code=attr['station'],
                            latitude=attr['latitude'],
                            longitude=attr['longitude'],
                            elevation=attr['elevation'],
                            creation_date=UTCDateTime(2016, 1, 2),
                            site=Site(name=attr['description']))
        
        attr['location'] = alphanum[locationindex]

        # accurate location if available
        if coordinates is not None:
            attr['latitude'] = latitudes[i]
            attr['longitude'] = longitudes[i]
            attr['elevation'] = elevations[i]

        # strain channel if strain file is provided
        if strainfile is not None:
            
            attr['npts'] = len(strain_array[:,i])
            attr['channel'] = '%s%s%s'%(band_code,strain_code,orientation_code) 
            stream += Trace(data=strain_array[:,i],header=attr)

            cha = Channel(
                code=attr['channel'], 
                location_code=attr['location'],
                latitude=attr['latitude'],
                longitude=attr['longitude'],
                elevation=attr['elevation'],
                depth=attr['depth'],
                azimuth=attr['azimuth'],
                dip=attr['dip'],
                sample_rate=attr['sampling_rate'])

            # Now tie it all together.
            cha.response = strain_response
            sta.channels.append(cha)

        # skip the velocity channel if velocity file is not provided
        if velocityfile is None:
            continue

        attr['npts'] = len(strain_array[:,i])
        attr['channel'] = '%sH%s'%(band_code,orientation_code)

        stream += Trace(data=velocity_array[:,i],header=attr)

        cha = Channel(
            code=attr['channel'], 
            location_code=attr['location'],
            latitude=attr['latitude'],
            longitude=attr['longitude'],
            elevation=attr['elevation'],
            depth=attr['depth'],
            azimuth=attr['azimuth'],
            dip=attr['dip'],
            sample_rate=attr['sampling_rate'])

        # Now tie it all together.
        cha.response = velocity_response
        sta.channels.append(cha)

        
    if sta is not None:
        net.stations.append(sta)
        
    inv.networks.append(net)

    # Write station inventory to a StationXML file. We also force a validation against
    # the StationXML schema to ensure it produces a valid StationXML file.
    inv.write(metadatastationxml, format="stationxml", validate=True)
    print(inv)

    # Write trace stream to miniseed file. We force record length and encoding style 
    # for consistency with SeisComP
    print(stream.__str__(extended=True))
    stream.write(datamseed, format="MSEED", reclen=512, encoding='FLOAT64')

    # And finally write SeisComP station configuration to an SCXML file.
    ## begin bindings
    if confxml is not None:
        binding ="""<?xml version="1.0" encoding="UTF-8"?>
        <seiscomp %s>
            <Config>
        """%seiscompschema
        
        ## indiv station bindings
        for station in stations:
            binding += bindingtemplate%(attr['network'],
                                        station,attr['network'],
                                        station,attr['network'],
                                        station,attr['network'],
                                        station,band_code,
                                        strain_code,orientation_code)
        ## finish bindings    
        binding += '''</Config>
        </seiscomp>'''

        with open(confxml, "w") as binding_file:
            binding_file.write(binding)
    
    return inv,stream,binding
