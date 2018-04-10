"""
Functions to calculate total score for the molecular candidate.

@author: Dr. Ivan Laponogov
"""


def total_multiplicative_score(scores):
    result=1.0;
    for s in scores:
        result=result*scores[s];
    return result;
        
def total_frag_score(scores):
    if 'Frag' in scores:
        return scores['Frag'];
    else:
        return 1.0;

def total_fpt_score(scores):
    if 'FPT' in scores:
        return scores['FPT'];
    else:
        return 1.0;

def total_fptfrag_score(scores):
    result=1.0;
    if 'FPT' in scores:
        result=result*scores['FPT'];
    if 'Frag' in scores:
        result=result*scores['Frag'];
    return result;


if __name__=='__main__':
    scores={'FPT':5,'Frag':3,'CFM':8};
    scores2={'FPT':5};
    print(total_multiplicative_score(scores));
    print(total_fpt_score(scores));
    print(total_fpt_score(scores2));
    print(total_fptfrag_score(scores));
    print(total_fptfrag_score(scores2));
