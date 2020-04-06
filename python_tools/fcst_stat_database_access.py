#  fcst_stat_database_access.py
#  Created: 04-2020    by:  Bryan Karpowicz (GMAO/GESTAR/USRA)
#  basically a clone similar to obs_databse_access.py Created: 01-2017    by:  Will McCarty (GMAO)
#
# This python module interfaces with the SEMPERPY/GMAOPY PSQL database 
# utilized on NCCS/Discover.  This module is not intended to plot the 
# forecaset stats; its aim is to directly interface/return values and statistics
# from the database so that the user can generate their own plots (or wrappers)
#
# Dependencies:  Numpy, copy, pprint, psycopg2 
#

import numpy as np
import psycopg2 as psql

class experiment():
# Class experiment 
#
#  Description:  Class (interface and subroutines) used to connect to and query
#    the SEMPERPY/GMAOPY database on discover. 
#  Notes: 
#    - Currently, the default is the user database, though it should work with 
#      the OPS database as long as it is explicitly specified in the 
#      initialization - optional flag database='ob_ops' (I think).
#  Inputs (upon initialization):
#    exp:       experiment ID used to populate the database
#    database:  database name (optional, ob_exp by default)
#    startdate: YYYYMMDDHH start date (optional, will determine first date for
#                  experiment otherwise)
#    enddate:   YYYYMMDDHH end date (optional, will determine last date for
#                  experiment otherwise)
#  Returns:
#    self instance: variable containing the conenction and used to query
#                    for initialized experiment
#  Usage:
#    import obs_database_access as fcst_db
#
#    experiment = fcst_db.experiment('bmk_ana.21z', startdate=2018080100, enddate=2018093000)
#
#    height_500mb_tau_24hrs_northern_hemisphere = experiment.get_stats(500, 24, 'n.hem','cor','h')
#
#    # will get you a numpy array of the 500mb height anomaly for 24 hour forcasts in the nothern hemisphere

    def __init__(self, exp, database='ob_exp', verify='gmao',startdate=None, enddate=None, verbose_query=None):

        self.experiment = exp
        print('Initializing as ',exp)
        self.database = database
        self.verify = verify 
        if (database == 'ob_exp'):
           self.connect_string = "dbname='semper' user='gmao_user' host='edb1' "
        elif (database == 'ob_ops'):
           self.connect_string = "dbname='gmao_stats' user='gmao_user' host='edb1' "
        else:
           print('Warning, not DB ob_exp or ob_ops, using dbname \'semper\'')
           self.connect_string = "dbname='semper' user='gmao_user' host='edb1' "

        self.con = psql.connect(self.connect_string)

        if (startdate and enddate):
            self.startdate = startdate
            self.enddate = enddate
        elif (startdate or enddate):
            print('warning, only startdate or enddate set, determining it from DB')
            self.startdate, self.enddate = self.get_min_max_date()
        else:
            self.startdate, self.enddate = self.get_min_max_date()
        print('experiment start, end dates:',self.startdate, self.enddate)
        self.verbose_query = verbose_query

    def get_pressure_levels(self):
# Function exp.get_pressure_levels()
#  Description: Get the pressure levels available from the experiment's database entry
#  Inputs: None
#  Returns: presure levels [hPa]
        
        cur = self.con.cursor()
        print("Determining levels from database")
        query = " SELECT distinct level from fc_exp.v_view where date = {} and step = 0 and expver = '{}' and verify = '{}' ORDER BY level ".format(self.startdate, self.experiment, self.verify)
        cur.execute(query)
        val = cur.fetchall()
        levs = np.array(val)[:,0]
        levs = levs[::-1]
        msk = (levs <= maxlev) & (levs >= minlev)
        levs = levs[msk]
        return( levs )
    def get_stats(self, lev, stp,  dom, statistic, variable):
        cur = self.con.cursor()
        query = " SELECT value from fc_exp.v_view where date >= {} and date <= {} and level = {} and step = {} and domain_name = '{}' and statistic = '{}' and expver = '{}' and variable = '{}' and verify = '{}' ORDER BY date ".format(self.startdate, self.enddate, lev, stp, dom, statistic, self.experiment,variable, self.verify)
        cur.execute(query)
        val = np.array(cur.fetchall())
        return(val)
    def set_daterange(self, startdate, enddate):
# Function exp.set_daterange()
#  Description:  Alter the start/end dates set upon initialization.
#  Notes:
#    - Should complain if only one date is set.  Perhaps in the future
#      make dates optional so that only one can be adjusted on the fly
#  Inputs: 
#    startdate: YYYYMMDDHH start date 
#    enddate:   YYYYMMDDHH end date 
#  Returns: 
#    Nothing
#  Example:
#    #continuing from initialization example above)
#    exp.set_daterange(2006070100,2006073118)
        self.startdate = startdate
        self.enddate = enddate
        print('experiment start, end dates set to:',self.startdate, self.enddate)

    def get_min_max_date(self):
# Function exp.get_min_max_date()()
#  Description:  Determines the earliest and latest date for the experiment
#    in the database.
#  Notes:
#    - Used internally to set start/end date if not specified upon 
#        initialization
#  Inputs: 
#    Nothing
#  Returns: 
#    startdate:  YYYYMMDDHH for earliest date
#    enddate:    YYYYMMDDHH for latest date
#  Example:
#    #continuing from initialization example above)
#    start_date, end_date = exp.get_min_max_date() 
	    
        cur = self.con.cursor()
        query = "SELECT MIN(date), MAX(date) from {}.v_view where  expver = \'{}\'".format(self.database, self.experiment)
        cur.execute(query)
        val = cur.fetchall()
        val = val[0]
        startdate, enddate = val
        return(startdate, enddate)
