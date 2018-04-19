MGF_DEFAULT_PARAMS={
                    'PEPMASS': '',
                    'MSLEVEL': '',
                    'CHARGE': ''
}

CD_DEFAULT_PARAMS={
                'db_molecule_name': '',
                'exactmass': '',
                'formula': '',
                'fpt_0': '',
                'fptcount': '',
                'global_index': '',
                'inchi': '',
                'level': '1',
                'mode': '',
                'peaks': '',
                'charge': '',
                'ion_type': ''
                }
# level: int, mode: int, collision_energy: float
# these parameters will be loaded by chemdistiller.msspectra.spectrum as specific types
#'mode': '',
#'collision_energy': '',

CD_DEFAULT_SUB_PARAMS={
                    'charge': '',
                    'collision_energy': '',
                    'collision_record': '',
                    'dbsource': '',
                    'exactmass': '',
                    'formula': '',
                    'inchi': '',
                    'level': '2',
                    'mode': '',                 # compulsory parameter for level 1
                    'precursor_ion': '',
                    'precursor_mz': '',
                    'peaks':''
}
# this chunk begins with 'spectrum' and ends with 'end'
# sub peaks begins with 'peaks' and end with 'end'

#this is a map of keywords
