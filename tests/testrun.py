import sys;

if sys.byteorder!='little':
    print('Only little endian machines currently supported! bye bye ....');
    quit();
    
import os;    

sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../chemdistiller")));


from chemdistiller.scorers.formula import FormulaScorer;
from chemdistiller.scorers.element import ElementScorer;
from chemdistiller.chemdb.manager import DBManager;
from chemdistiller.filters.formula import FormulasFilter;

with DBManager() as DataBaseManager:
     test_formulas_filter=FormulasFilter('H4C4SO,H4C3SiO2');
     test_formula_scorer=FormulaScorer();
     test_formula_scorer.configure_scorer('H4C4SO,H4C3SiO2',[0.8,0.5],1.0);
         
     test_elements_scorer=ElementScorer();
     test_elements_scorer.configure_scorer({'S':0.1,'O':1.0,'Si':0.1});
     print(DataBaseManager.db_indexes_from_db_names(['PubChem','EMolecules']));
         
     result=DataBaseManager.query_by_mz_scored(100, 20, charge=0, db_indexes=DataBaseManager.db_indexes_from_db_names(['EMolecules','PubChem'],case_sensitive=True), filters=[test_formulas_filter], scorers=[test_formula_scorer,test_elements_scorer], required_fields=set(['Formula']), results_limit=20, save_memory=False);
     
     print(len(result.mol_list));
     for i in result.mol_list:
         print(i['MZ'],i['DBIndex'],i['Formula'],DataBaseManager.data_bases[i['DBIndex']].db_name,i['Scores']);
     
