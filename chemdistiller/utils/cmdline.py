# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 12:43:29 2017

@author: Dr. Ivan Laponogov
"""
import sys;
import os;
import json;
import numpy as np;

class AttributeCollection:
    pass;
    
class AllValues:
    pass;
    
class UnrecognisedArguments(Exception):
    pass

class WrongArgumentDefinition(Exception):
    pass

class WrongConfigurationDefinition(Exception):
    pass

def de_format(string):
    string=string.replace('\n', ' ').replace('\t',' ');
    while '  ' in string:
        string=string.replace('  ',' ');
    return string.strip();

def format_to_width(string, width, add_offset):
    string=string.split();
    out_strings=[];
    current_string='';
    for s in string:
        if len(current_string)+len(s)+1>width:
            out_strings.append(current_string);
            current_string='%s%s'%((' '*add_offset), s);
        else:
            if current_string!='':
                current_string=" ".join((current_string,s));
            else:
                current_string=s;
    if current_string!='':
        out_strings.append(current_string);
    return out_strings;

class Option:
    
    @staticmethod
    def _equallists(targets_a, targets_b):
        if len(targets_a)!=len(targets_b):
            return False;
        ta=sorted(targets_a);
        tb=sorted(targets_b);
        for i in range(len(ta)):
            if ta[i]!=tb[i]:
                return False;
        return True;
    
    def _check_targets(self, targets):
        option=self._option.lstrip('-');
        if option in targets:
            if not self._equallists(self._targets, targets[option]):
                raise WrongConfigurationDefinition('Error! Option "%s" has different targets defined in different places! This will cause the confusion for interpreter!'%self._option);
        else:
            targets[option]=self._targets;
        if self._values:
            for value in self._values:
                if isinstance(value, Value):
                    if value._parameters:
                        for parameter in value._parameters:
                            parameter._check_targets(targets);
        
        
    def __init__(self, option, help='', values=[None], is_list=False, type=None, conditions=[], targets=[], optional=True):
        self._option='';
        if ',' in option:
            option=option.split(',');
            for suboption in option:
                if suboption.startswith('--'):
                    self._option=suboption;
                elif suboption.startswith('-'):
                    self._short_option=suboption;
                else:
                    raise WrongArgumentDefinition('Multiple names for option are not allowed! Option: "%s"'%option);
        else:
            self._option=option;    
            self._short_option='';

        if self._option=='':
            raise WrongArgumentDefinition('Empty full name for option is not allowed! Option: "%s"'%option);
        
        self._is_list=is_list;
        self._help=help;
        self._values=values;
        self._type=type;
        self._conditions=conditions;
        self._targets=targets;
        self._optional=optional;
        #if values:
        #    self.current_value=values[0];

    def _value_as_text(self, value):
            if isinstance(value, Value):
                return self._value_as_text(value._value);
            if isinstance(value,str):
                return "'%s'"%value;
            else:
                return '%s'%value;
    
    def _get_values_list(self):
        lst=[];
        for value in self._values:
            lst.append(self._value_as_text(value));
        return lst;
        
    def _get_arg_values_def(self):
        if None in self._values:
            return self._option.lstrip('-').upper();
        else:
            return '{%s}'%('|'.join(self._get_values_list()));
    
    def _list_arguments(self, args):
        args.append(self._option);
        if self._values:
            for value in self._values:
                if isinstance(value, Value):
                    if value._parameters:
                        for parameter in value._parameters:
                            parameter._list_arguments(args);
    
    def _collect_arguments(self, args):
        if self._short_option!='':
            argname=', '.join((self._short_option, self._option));
        else:
            argname=self._option;
        
        if self._values:
            if self._option.startswith('--'):
                argname='%s %s'%(argname, self._get_arg_values_def());
            
        if self._optional:
            argname='[%s]'%argname;
        
        if not argname in args:
            args.append(argname);

        if self._values:
            for value in self._values:
                if isinstance(value, Value):
                    value._collect_arguments(args);
                    
    def _expand_targets(self, deftarget):
        self.children={};
        if self._values:
            for value in self._values:
                if isinstance(value, Value):
                    self.children[value._value]=value;
                    value.parent=self;
                    value._expand_targets('%s/%s'%(deftarget,self._option.lstrip('-')))

    def _get_description(self, include_targets, line_width, offset):
        
        helpstr=self._help;
        helpstr='%s  : %s'%(str(self), helpstr)
        if self._values:
            if not (self._values[0] is None):
                helpstr+=' Default value: %s.\n'%self._value_as_text(self._values[0]);
        
        helpstr=format_to_width(de_format(helpstr), line_width-offset-4, 4);
        description=(' '*(offset+4))+('\n%s'%(' '*(offset+4))).join(helpstr);
            
        description+='\n';
        
        if self._values:
            for value in self._values:
                if isinstance(value, Value):
                    description+='%s\n'%value._get_description(include_targets, line_width, offset+4);
        #print(include_targets);
        if include_targets:
            description+='%sTargets: %s\n'%(' '*(offset+3),self._targets);
        
        return description;
        
    def __repr__(self):
        if self._short_option!='':
            option='%s, %s'%(self._short_option, self._option);
        else:
            option=self._option;
            
        return '%s %s'%(option, self._get_arg_values_def());
    
    
    def _process_argvalue(self, argvalue):
        if isinstance(argvalue, list):
            for i in range(len(argvalue)):
                argvalue[i]=self._process_argvalue(argvalue[i]);
        else:
            if not (self._type is None):
                try:
                    argvalue=self._type(argvalue);
                except:
                    raise WrongArgumentDefinition('Wrong type! Argument "%s" requires %s!'%(self._option, self._type));
            if self._values:
                if not (None in self._values):
                    accepted=False;
                    for value in self._values:
                        if isinstance(value, Value):
                            value=value._value;
                        if value==argvalue:
                            accepted=True;
                            break;
                    if not accepted:
                        raise WrongArgumentDefinition('Wrong value! Argument "%s" does not accept value %s!'%(self._option, argvalue));

            if self._conditions:            
                self.__dict__[self._option.lstrip('--')]=argvalue;
                for condition in self._conditions:
                    if not condition[1](self):
                        raise WrongArgumentDefinition('Wrong value! For argument "%s" condition %s is not satisfied!'%(self._option, condition[0]));
        
        return argvalue;
    
    
    def _process_subvalues(self, argvalue, parsed_args, parameters, param_path, root):
        if self._values:
            for value in self._values:
                 if isinstance(value, Value):
                      if argvalue==value._value:
                          value._parse_args(parsed_args, parameters, param_path, root);
                          
    
    def _parse_args(self, parsed_args, parameters, param_path, root):
        already_set=root._already_set(self._targets, self._option.lstrip('-'));
        if already_set[0]:
            argvalue=already_set[1];
            
        else:
            if self._option.startswith('--'):
                
                if self._short_option!='':
                    while self._short_option in parsed_args:
                        parsed_args[parsed_args.index(self._short_option)]=self._option;
                
                argvalue=None;
    
                if not (self._option in parsed_args):
                    if self._optional:
                        if self._values:
                            argvalue=self._values[0];
                            if isinstance(argvalue, AllValues):
                                argvalue=self._values[1:];
    
                            if isinstance(argvalue, list):
                                for i in reversed(range(len(argvalue))):
                                    if argvalue[i] is None:
                                        del argvalue[i];
                                    else:    
                                        if isinstance(argvalue[i], Value):
                                            argvalue[i]=argvalue[i]._value;
                                        if isinstance(argvalue[i], str):
                                            if argvalue[i].startswith('"') and argvalue[i].endswith('"'):
                                                argvalue[i]=argvalue[i][1:-1];
                                            elif argvalue[i].startswith("'") and argvalue[i].endswith("'"):
                                                argvalue[i]=argvalue[i][1:-1];
                            else:        
                                if isinstance(argvalue, Value):
                                    argvalue=argvalue._value;
                                if isinstance(argvalue, str):
                                    if argvalue.startswith('"') and argvalue.endswith('"'):
                                        argvalue=argvalue[1:-1];
                                    elif argvalue.startswith("'") and argvalue.endswith("'"):
                                        argvalue=argvalue[1:-1];
                    else:
                        if not root._initializing:
                            raise WrongArgumentDefinition('Argument "%s" is not optional and must be provided!'%self._option)
                else:
                    while (self._option in parsed_args):
                        index=parsed_args.index(self._option);
                        if self._values:
                            try:
                                argvalue=parsed_args[index+1];
                                del parsed_args[index+1];
                                del parsed_args[index];
                                if self._is_list:
                                    argvalue=argvalue.split(',');
                                    for i in range(len(argvalue)):
                                        if argvalue[i].startswith('--'):
                                            raise WrongArgumentDefinition('Argument "%s" requires a value!'%self._option);
                                        if argvalue[i].startswith('-'):
                                            f=float(argvalue[i]); #Test if negative number
                                else:
                                    if argvalue.startswith('--'):
                                        raise WrongArgumentDefinition('Argument "%s" requires a value!'%self._option);
                                    if argvalue.startswith('-'):
                                        f=float(argvalue); #Test if negative number
                            except:
                                raise WrongArgumentDefinition('Argument "%s" requires a value!'%self._option);
    
                            if isinstance(argvalue, list):
                                for i in range(len(argvalue)):
                                    if isinstance(argvalue[i], str):
                                        if argvalue[i].startswith('"') and argvalue[i].endswith('"'):
                                            argvalue[i]=argvalue[i][1:-1];
                                        elif argvalue[i].startswith("'") and argvalue[i].endswith("'"):
                                            argvalue[i]=argvalue[i][1:-1];
                            else:
                                if isinstance(argvalue, str):
                                    if argvalue.startswith('"') and argvalue.endswith('"'):
                                        argvalue=argvalue[1:-1];
                                    elif argvalue.startswith("'") and argvalue.endswith("'"):
                                        argvalue=argvalue[1:-1];
    
                            argvalue=self._process_argvalue(argvalue);
                
            else:
                argvalue=None;
                if len(parsed_args)>0:
                            try:
                                argvalue=parsed_args[0];
                                del parsed_args[0];
                                if self._is_list:
                                    argvalue=argvalue.split(',');
                                    for i in range(len(argvalue)):
                                        if argvalue[i].startswith('--'):
                                            raise WrongArgumentDefinition('Argument "%s" requires a value!'%self._option);
                                        if argvalue[i].startswith('-'):
                                            f=float(argvalue[i]); #Test if negative number
                                else:
                                    if argvalue.startswith('--'):
                                        raise WrongArgumentDefinition('Argument "%s" requires a value!'%self._option);
                                    if argvalue.startswith('-'):
                                        f=float(argvalue); #Test if negative number
                            except:
                                raise WrongArgumentDefinition('Argument "%s" requires a value!'%self._option);
    
                            if isinstance(argvalue, list):
                                for i in range(len(argvalue)):
                                    if isinstance(argvalue[i], str):
                                        if argvalue[i].startswith('"') and argvalue[i].endswith('"'):
                                            argvalue[i]=argvalue[i][1:-1];
                                        elif argvalue[i].startswith("'") and argvalue[i].endswith("'"):
                                            argvalue[i]=argvalue[i][1:-1];
                            else:
                                if isinstance(argvalue, str):
                                    if argvalue.startswith('"') and argvalue.endswith('"'):
                                        argvalue=argvalue[1:-1];
                                    elif argvalue.startswith("'") and argvalue.endswith("'"):
                                        argvalue=argvalue[1:-1];
    
                            argvalue=self._process_argvalue(argvalue);
                
                else:
                    if self._optional  or root._initializing:
                        if self._values:
                            argvalue=self._values[0];
                            if isinstance(argvalue, AllValues):
                                argvalue=self._values[1:];
    
                            if isinstance(argvalue, list):
                                for i in reversed(range(len(argvalue))):
                                    if argvalue[i] is None:
                                        del argvalue[i];
                                    else:    
                                        if isinstance(argvalue[i], Value):
                                            argvalue[i]=argvalue[i]._value;
                                        if isinstance(argvalue[i], str):
                                            if argvalue[i].startswith('"') and argvalue[i].endswith('"'):
                                                argvalue[i]=argvalue[i][1:-1];
                                            elif argvalue[i].startswith("'") and argvalue[i].endswith("'"):
                                                argvalue[i]=argvalue[i][1:-1];
                            else:        
                                if isinstance(argvalue, Value):
                                    argvalue=argvalue._value;
                                if isinstance(argvalue, str):
                                    if argvalue.startswith('"') and argvalue.endswith('"'):
                                        argvalue=argvalue[1:-1];
                                    elif argvalue.startswith("'") and argvalue.endswith("'"):
                                        argvalue=argvalue[1:-1];
                        
                    else:
                        if not root._initializing:
                            raise WrongArgumentDefinition('Argument "%s" is not optional and must be provided!'%self._option) 
        self._current_value=argvalue;
        if not already_set[0]:
            root._set_arg_value(self._targets, self._option.lstrip('-'), argvalue);
        if isinstance(argvalue, list):
            for value in argvalue:
                self._process_subvalues(value, parsed_args, parameters, param_path, root);
        else:
            self._process_subvalues(argvalue, parsed_args, parameters, param_path, root);
        
class Value:
    def __init__(self, value, help='', parameters=[]):
        self._value=value;
        self._help=help;
        self._parameters=parameters;
    
    def _collect_arguments(self, args):    
        for parameter in self._parameters:
            parameter._collect_arguments(args);


    def _expand_targets(self, subfolder):
        self.children={};
        for parameter in self._parameters:
            self.children[parameter._option.lstrip('-')]=parameter;
            parameter.parent=self;
            if not parameter._targets:
                parameter._targets=[subfolder+'_'+self._value];
            parameter._expand_targets(subfolder+'_'+self._value);
                
    def _get_description(self, include_targets, line_width, offset):

        helpstr=de_format(self._help)

        if self._parameters:
            helpstr='For %s (%s) option(s):'%(str(self), helpstr);
        else:
            helpstr='%s (%s)'%(str(self), helpstr);
            
        helpstr=format_to_width(helpstr, line_width-offset-4,4);
        
        description='\n'+(' '*(offset+4))+('\n%s'%(' '*(offset+4))).join(helpstr)+'\n';
        
        if self._parameters:
            for parameter in self._parameters:
                description+='\n%s'%parameter._get_description(include_targets, line_width, offset+4);
        
        return description;

    def __repr__(self):
        if isinstance(self._value,str):
            return "'%s'"%self._value;
        else:
            return '%s'%self._value;    
        
    def _parse_args(self, parsed_args, parameters, param_path, root):
        for parameter in self._parameters:
            parameter._parse_args(parsed_args, parameters, param_path, root);
        
#-------------------------------------------                
    
class OptionsHolder:


    def _resolve_path(self, path):
        path=path.lstrip('/');
        result=self.parameters;
        if path!='':
            path=path.split('/');
        for subpath in path:
            if not subpath in result:
                result[subpath]={};
            result=result[subpath];
        return result;
    
    def _already_set(self, targets, option):
        result=False;
        value=None;
        for target in targets:
            subpath=self._resolve_path(target);
            if option in subpath:
                result=True;
                value=subpath[option];
                break
        return (result, value);
    
    def _set_arg_value(self, targets, option, argvalue):
        for target in targets:
            self._resolve_path(target)[option]=argvalue;
            
    def _generate_arguments_description(self, include_targets):
        description='Positional arguments:\n\n';
        
        for option in self._configuration:
            if not option._option.startswith('--'):
                description+='%s\n'%option._get_description(include_targets, self._default_line_width, 0);
                
        description+='\nOptional arguments:\n\n   -h, --help  : Show this help message and exit\n\n';
        
        for option in self._configuration:
            if option._option.startswith('--'):
                description+='%s\n'%option._get_description(include_targets, self._default_line_width, 0);
        
        return description;
           


    #------------------Finished--------------------------

    def __init__(self, docstring, configuration):
        program_description=[];
        description_epilog=['=============================================================================='];
        docstring=docstring.split('\n');
        stat=0;
        for s in docstring:
            if s.startswith('run python'):
                stat=1;
            elif s.startswith('Copyright:'):
                stat=2;
            #else:
            #    stat=3;
            if stat==0:
                if '***********' in s:
                    s='==============================================================================';
                program_description.append(s);
            elif stat==2:
                description_epilog.append(s);
            
        self._default_line_width=79;        
        self.program_description='\n'.join(program_description);
        self.description_epilog='\n'.join(description_epilog);
        self._configuration=configuration;
        self._expand_targets();
        self.exec_name='';
        self._check_sanity();
        self._initializing=True;
        self._process_parsed_args();
        self._initializing=False;
        

    
    def _process_parsed_args(self, parsed_args=None):
        if parsed_args is None:
            parsed_args=[];
        self.parameters={};
        for option in self._configuration:
            option._parse_args(parsed_args, self.parameters, '', self);
        if parsed_args:
            raise UnrecognisedArguments('Arguments not recognised: %s'%parsed_args)

    def _generate_arguments_list(self):
        offset=' '*(len(self.exec_name)+1);
        args=[];
        for option in self._configuration:
            option._collect_arguments(args);
        if args:
            args.insert(0,'[-h, --help]');
            arguments_list=(',\n%s'%offset).join(args);
        return arguments_list;

    def parse_command_line_args(self):
        self._exec_name=sys.argv[0];
        args=sys.argv[1:];

        if not args: #print usage if no arguments supplied and exit
            print(self.format_usage());
            print('\n');
            print(self.description_epilog);
            sys.exit(0);
        elif '-h' in args or '--help' in args: #print help and exit if help argument supplied
            print(self.format_help());
            sys.exit(0);
        else:
            self.parse_args(args);
       
    def _expand_targets(self):
        self.children={};
        for option in self._configuration:
            option.parent=self;
            self.children[option._option.lstrip('-')]=option;
            if not option._targets:
                option._targets=['/'];
            option._expand_targets('');    
       
    def parse_args(self, args=[]):
        #print(args)
        if isinstance(args, str):
            args=args.split();
        if args:
            for i in reversed(range(len(args)-1)):
                if args[i].endswith(','):
                    args[i]=args[i]+args.pop(i+1);
        self._process_parsed_args(args);
    
    def format_help(self, include_targets=False):
        help_string="\n\n".join([self.format_usage(), self._generate_arguments_description(include_targets), self.description_epilog]);
        return help_string;
    
    def format_usage(self):
        return 'Usage:\n\n%s %s'%(os.path.basename(self._exec_name), self._generate_arguments_list());
    
    def _format_parameter_settings(self, parameters, prefix=''):
        result=[];
        for key in sorted(parameters.keys()):
            if isinstance(parameters[key], dict):
                if prefix!='':
                    sub='.'.join(prefix,key);
                else:
                    sub=key;
                result.append(self._format_parameter_settings(parameters[key], sub));
            else:
                keyvalue=parameters[key];
                if isinstance(keyvalue, str):
                    keyvalue='"%s"'%keyvalue;
                if prefix=='':
                    result.append('%s = %s'%(key, keyvalue));
                else:
                    result.append('%s.%s = %s'%(prefix, key, keyvalue));
                        
                    
        return '\n'.join(result);
    
    def format_parameters(self):
        return '\n'.join(('Current parameters:', self._format_parameter_settings(self.parameters,'')));
            
        
    def export_as_json(self):
        return json.dumps(self.parameters, indent=4, separators=(',', ': '));

    def _check_sanity(self):
        targets={};
        for option in self._configuration:
            option._check_targets(targets);

            
    #def store_to_hdf5(self, h5file):
    #    pass
    
    #def restore_from_hdf5(self, h5file):
    #    pass
    def _list_arguments(self):
        args=[];
        for option in self._configuration:
            option._list_arguments(args);
            
        unique_args=[];
        for arg in args:
            if not arg in unique_args:
                unique_args.append(arg);
        return unique_args;
    

    def _de_tree(self, args, dictionary):
        for key in dictionary.keys():
            if isinstance(dictionary[key], dict):
                self._de_tree(args, dictionary[key]);
            else:
                args[key]=dictionary[key];
            
    
    def import_from_dict(self, dictionary):
        lst=self._list_arguments();
        args={};
        self._de_tree(args, dictionary);

        argarray=[];

        for option in lst:
            key=option.lstrip('-');
            if key in args:
                if option.startswith('--'):
                    argarray.append(option);
                argarray.append(str(args[key]));
        
        self.parse_args(argarray);
    
        
    def import_from_json(self, js=''):
        self.import_from_dict(json.loads(js));
        
        
        
if __name__=='__main__':
    
    if sys.byteorder!='little':
        print('Only little endian machines currently supported! bye bye ....');
        quit();

    module_path = os.path.abspath('%s/..'%os.path.dirname(os.path.realpath(__file__)));
    sys.path.append(module_path);

    import procconfig as pf;
    optionsholder=OptionsHolder(__doc__, pf.PeakAlign_options);
    
    optionsholder._exec_name=os.path.basename(__file__);
    
    args='--peakalignmethod NN test.hdf5 out.hdf5 --maxpeakshift 0.3 --units Da'    
    
    optionsholder.parse_args(args);
    
    print(optionsholder.format_parameters());
    

