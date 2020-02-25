#!/usr/bin/env python3

import sys
import re
import numpy as np
import datetime

# Version tag (here comes the version tag)
#
# (here comes the COMMIT info)
# (here comes the DATE info)
try:
    version_tag
except:
    class version_tag:
        COMMIT="Untracked"
        DATE="No date"

# Define classes for atom, bond...
class atom:
    def __init__(self,*args):
#    def __init__(self,iat,attype,ires,resname,atname,chgr,q,mass):
        if len(args) == 8:
            self.iat     = int(args[0])
            self.attype  = args[1]
            self.ires    = int(args[2])
            self.resname = args[3]
            self.atname  = args[4]
            self.chrg    = int(args[5])
            self.q       = float(args[6])
            self.mass    = float(args[7])
            self.V       = 0.0  # sigma(cr=2,3) or C6(cr=1)
            self.W       = 0.0  # epsilon(cr=2,3) or C12(cr=1)
        elif len(args) == 7:
            # E.g. IB+ ion
            self.iat     = int(args[0])
            self.attype  = args[1]
            self.ires    = int(args[2])
            self.resname = args[3]
            self.atname  = args[4]
            self.chrg    = int(args[5])
            self.q       = float(args[6])
            self.mass    = 0.0  # To be set..
            self.V       = 0.0  # sigma(cr=2,3) or C6(cr=1)
            self.W       = 0.0  # epsilon(cr=2,3) or C12(cr=1)
        else:
            raise TypeError('Wrong number of elements to instaciate atom. Expected 7 or 8, got '+str(len(args)))
    def setLJ(self,V,W):
        self.V       = float(V)
        self.W       = float(W)
        
class bond:
    def __init__(self):
        self.i1   = 0
        self.i2   = 0
        self.ft   = 0
        self.r0   = 0.0
        self.kb   = 0.0
        self.prms = ''        
    def setbond(self,i1,i2,ft,params):
        self.i1   = int(i1)
        self.i2   = int(i2)
        self.ft   = int(ft)
        self.prms = params
    def entryline(self):
        return "%5i %5i %2i %s"%(self.i1,self.i2,self.ft,self.prms)


class pair:
    def __init__(self):
        self.i1   = 0
        self.i2   = 0
        self.ft   = 0
        self.prms = ''        
    def setpair(self,i1,i2,ft,params):
        self.i1   = int(i1)
        self.i2   = int(i2)
        self.ft   = int(ft)
        self.prms = params
    def entryline(self):
        return "%5i %5i %2i %s"%(self.i1,self.i2,self.ft,self.prms)
        
class angle:
    def __init__(self):
        self.i1 = 0
        self.i2 = 0
        self.i3 = 0
        self.ft = 0
        self.prms = ''
    def setangle(self,i1,i2,i3,ft,params):
        self.i1   = int(i1)
        self.i2   = int(i2)
        self.i3   = int(i3)
        self.ft   = int(ft)
        self.prms = params
    def entryline(self):
        return "%5i %5i %5i %2i %s"%(self.i1,self.i2,self.i3,self.ft,self.prms)
        
class dihed:
    def __init__(self):
        self.i1 = 0
        self.i2 = 0
        self.i3 = 0
        self.i4 = 0
        self.ft = 0
        self.prms = ''
    def setdihed(self,i1,i2,i3,i4,ft,params):
        self.i1   = int(i1)
        self.i2   = int(i2)
        self.i3   = int(i3)
        self.i4   = int(i4)
        self.ft   = int(ft)
        self.prms = params
    def entryline(self):
        return "%5i %5i %5i %5i %2i %s"%(self.i1,self.i2,self.i3,self.i4,self.ft,self.prms)

# INPUT PARSER
def get_args():
    
    # Options and their defaults 
    final_arguments = dict()
    final_arguments["-f"]="input.dat"
    final_arguments["-s"]="none"
    final_arguments["-h"]=False
    # Description of the options
    arg_description = dict()
    arg_description["-f"] ="Name of the input topology (processed)"
    arg_description["-s"] ="Swapp file (if not equal none)"
    arg_description["-h"] ="Show this help"
    # Type for arguments
    arg_type = dict()
    arg_type["-f"] ="char"
    arg_type["-s"] ="char"
    arg_type["-h"]    ="-"
    
    # Get list of input args
    input_args_list = []
    iarg = -1
    for s in sys.argv[1:]:
        # get -flag [val] arguments 
        if s[0]=="-":
            iarg=iarg+1
            input_args_list.append([s])
        else:
            input_args_list[iarg].append(s)
            
    # Transform into dict. Associtaing lonely flats to boolean   
    input_args_dict=dict()
    for input_arg in input_args_list:
        if len(input_arg) == 1:
            # Boolean option. Can be -Bool or -noBool
            input_arg.append(True)
            if input_arg[0][1:3] == "no":
                input_arg[0] = "-" + input_arg[0][3:]
                input_arg[1] = not input_arg[1]
        elif len(input_arg) != 2:
            raise BaseException("Sintax error. Too many arguments")

        input_args_dict[input_arg[0]] = input_arg[1]
    
    for key,value in input_args_dict.items():
        # Check it is allowed
        isValid = final_arguments.get(key,None)
        if isValid is None:
            raise BaseException("Sintax error. Unknown label: " + key)
        # If valid, update final argument
        final_arguments[key]=value
        
    if final_arguments.get("-h"):
        
        print("""
 ----------------------------------------
        procesed2standalone.py

   Convert topology based on [ Xtypes ]
   databases to explicitly defined 
   potential terms in the topology
   
   Version(GIT HASH): %s
   Date             : %s
 ----------------------------------------
        """%(version_tag.COMMIT,version_tag.DATE))
        print("    Options:")
        print("    --------")
        print('      {0:<10}  {1:^4}  {2:<41}  {3:<7}'.format("Flag","Type","Description","Value"))
        print('      {0:-<10}  {1:-^4}  {2:-<41}  {3:-<7}'.format("","","",""))
        for key,value in final_arguments.items():
            descr = arg_description[key]
            atype = arg_type[key]
            #atype=str(type(value)).replace("<type '","").replace("'>","")
            print('      {0:<10}  {1:^4}  {2:<41}  {3:<7}'.format(key, atype, descr, str(value)))
        print("")
        
        sys.exit()
        
    return final_arguments


#############################################################
#
#  MAIN
#
#############################################################
if __name__ == '__main__':

    # Get arguments
    args=get_args()
    topfile = args.get('-f')

    # Read input file
    with open(topfile) as f:
        topentry=[]
        for line in f:
            line=line.replace('\n','')
            line=line.replace('\t','  ')
            line=line.split(';')[0]
            if len(line.lstrip()) != 0:
                topentry.append(line)
                
    # Read swap info (TBD)
    swapfile = args.get('-s')
    if swapfile != 'none':
        with open(swapfile) as f:
            swap=dict()
            for line in f:
                line=line.replace('\n','')
                line=line.replace('\t','  ')
                line=line.split(';')[0]
                if len(line.lstrip()) != 0:
                    line=line.split('--')[0]
                    i,j = line.split()
                    swap[int(i)] = int(j)
    
    # Prepare arrays for molecule items
    atoms=[]
    bonds=[]
    pairs=[]
    angles=[]
    diheds=[]
    # Prepare dictionaries for data types database
    atom_prms=dict()
    bond_prms=dict()
    pair_prms=dict()
    angle_prms=dict()
    dihed_prms=dict()
    
    
    print('; Topology Generated with '+sys.argv[0])
    print('; Launched on: ',datetime.datetime.now())
    print('; Original topology file: '+topfile)
    print(';')
    
    # LOOP0: get molecules
    # Initialize
    molecules=[]
    section=''
    for line in topentry:
        # Use re to locate sections
        # We use a named group to individuate the section. This is done with
        # "?P<sect>" (sections are the blocks in parenthesis)
        pattern = r'(^( )*\[( )*(?P<section>[a-zA-Z_]+)( )*\]( )*$)'
        match = re.match(pattern,line)
        if match:
            section=match.group('section')
        
        elif section == 'molecules':
            line_strip = line.split(';')[0]
            if (line_strip.replace(' ','')) == 0:
                continue
            else:
                data=line_strip.split()
                molecules.append(data[0])
    
    # LOOP1: get atoms
    # Initialize
    section=''
    for line in topentry:
        # Use re to locate sections
        # We use a named group to individuate the section. This is done with
        # "?P<sect>" (sections are the blocks in parenthesis)
        pattern = r'(^( )*\[( )*(?P<section>[a-zA-Z_]+)( )*\]( )*$)'
        match = re.match(pattern,line)
        if match:
            section=match.group('section')
            continue
        
        elif section == 'atoms' and molecule_exists:
            # Get all items for the atom in data array
            data = line.split()
            # Form the atoms array with a new atom (class)
            # *data expands the array into the elements
            atoms.append(atom(*data))
            # Get LJ parameters from database...

        if section == 'moleculetype':
            molname = line.split()[0]
            if molname in molecules:
                molecule_exists=True
            else:
                molecule_exists=False
   

    # LOOP2: get dict types and used atomtypes
    # Initialize
    f = open('ff_gau.prm','w')
    print('! Non-bonded master equation for AmberFF;',file=f)
    print('NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.2',file=f)
    gau_bonds=[]
    gau_angles=[]
    gau_diheds=[]
    gau_diheds_amb=dict()
    section=''
    iat=0
    nshift=0
    molecule_exists=True
    for line in topentry:
        # Use re to locate sections
        # We use a named group to individuate the section. This is done with
        # "?P<sect>" (sections are the blocks in parenthesis)
        pattern = r'(^( )*\[( )*(?P<section>[a-zA-Z_]+)( )*\]( )*$)'
        match = re.match(pattern,line)
        if match:
            section=match.group('section')
            if section == 'atomtypes' or section == 'system' or section == 'molecules':
                print('\n%s'%(line))
            elif section == 'implicit_genborn_params':
                pass
            elif not 'type' in section and molecule_exists:
                print('\n%s'%(line))
            continue
            
        # Sections with parameters databases (types)
        #-------------------------------------------
        if section == 'defaults':
            data = line.split()
            nbfunc  = int(data[0])
            combrule= int(data[1])
            genpairs= data[2]
            fudgeLJ = float(data[3])
            fudgeQQ = float(data[4])
            print(line)
        elif section == 'atomtypes':
            # Get data
            data = line.split()
            # Select format based on the position of the ptype
            ptypes=['A','S','V','D']
            if data[4] in ptypes:
                # Most general (CHARMM,AMBER...)
                attype,atnum,atmass,atq,ptype,sigma,epsilon=data
                attypeb=attype
            if data[5] in ptypes:
                # OPLS: includes an attype to select the bonds (attypeb)
                attype,attypeb,atnum,atmass,atq,ptype,sigma,epsilon=data
            elif data[3] in ptypes:
                # Joyce: do not use atnum
                attype,atmass,atq,ptype,sigma,epsilon=data
                attypeb=attype
            # Reconstruct entry with Picky/Joyce expected format
            atom_prms[attype] = '  '.join([atmass,atq,ptype,sigma,epsilon,attypeb])
            # Only print used atoms
            if attype in [ atoms[i].attype for i in range(len(atoms)) ]:
                print('%-10s %10s %10s %1s %10s %10s'%(attype,atmass,atq,ptype,sigma,epsilon))
                #----------------------
                #Print gaussian params
                #----------------------
                rminh = float(sigma) * 2.**(1./6.) * 5.
                epsG  = float(epsilon)/4.184
                print('VDW %s  %10.3f  %10.3f'%(attype,rminh,epsG),file=f)
                #
                
                
        elif section == 'bondtypes':
            # Get data
            data = line.split()
            bondtype = '-'.join(data[:2])
            params   = '   '.join(data[3:])
            # Form dictionaries with data base
            # Entries are float (not list)
            # Generate new dict layer if still not created
            ft=int(data[2])
            if not ft in bond_prms:
                bond_prms[ft] = dict()
            # First check if already in to db
            if bondtype in bond_prms[ft]:
                print('; WARNING: multiple bond definition')
                print('; '+bondtype)
                print('; New parameters   : ',params)
                print('; Overridden params: ',bond_prms[ft][bondtype])
                print(';')
            else:
                bond_prms[ft][bondtype] = params
                
        elif section == 'pairtypes':
            # Get data
            data = line.split()
            pairtype = '-'.join(data[:2])
            params   = '   '.join(data[3:])
            # Form dictionaries with data base
            # Entries are float (not list)
            # Generate new dict layer if still not created
            ft=int(data[2])
            if not ft in pair_prms:
                pair_prms[ft] = dict()
            # First check if already in to db
            if pairtype in pair_prms[ft]: 
                print('; WARNING: multiple pair definition')
                print('; '+pairtype)
                print('; New parameters   : ',params)
                print('; Overridden params: ',pair_prms[ft][pairtype])
                print(';')
            else:
                pair_prms[ft][pairtype] = params
                
        elif section == 'angletypes':
            # Get data
            data = line.split()
            angletype = '-'.join(data[:3])
            params   = '   '.join(data[4:])
            # Form nested dictionaries with data base 
            # Entries are float (not list)
            # Generate new dict layer if still not created
            ft=int(data[3])
            if not ft in angle_prms:
                angle_prms[ft] = dict()
            # First check if already in to db
            if angletype in angle_prms[ft]:
                print('; WARNING: multiple angle definition')
                print('; '+angletype)
                print('; New parameters   : ',params)
                print('; Overridden params: ',angle_prms[ft][angletype])
                print(';')
            else:
                angle_prms[ft][angletype] = params
                
            angle_prms[angletype] = params
        elif section == 'dihedraltypes':
            # Get data
            data = line.split()
            # Ignore "exotic" dihedrals (see oplsaa.itp from VirtualChemistry: 
            #  [ dihedraltypes ]
            #  ;  i    j    k    l   func     coefficients
            #  ; Added DvdS for Quartz simulations
            #   SI   OS    1     0.000       3.766      3
            if data[2] in ['1','2','3','4']:
                continue
            dihedtype = '-'.join(data[:4])
            params   = '   '.join(data[5:])
            # Form dictionaries with data base
            # We need to make a nested dict:
            # dihed_prms[ft][dihedtype]
            # to differenciate entries for 
            # the same quarted for different ft
            # (e.g. proper and improper)
            # This is more general, but 
            # NOTE: for ft=9, there can be several
            #       dihedral parameters associated
            #       to one quartet
            # So, we use lists as entries for the dict
            
            # Generate new dict layer if still not created
            ft=int(data[4])
            if not ft in dihed_prms:
                dihed_prms[ft] = dict()
            
            # First check if already in to db
            if dihedtype in dihed_prms[ft]:
                if ft != 9:
                    print('; WARNING: multiple dihedral definition for type: '+str(ft))
                    print('; '+dihedtype)
                    print('; New parameters   : ',params)
                    print('; Overridden params: ',dihed_prms[ft][dihedtype][0])
                    print(';')
                # Add new params
                dihed_prms[ft][dihedtype].append(params)
                
            else:
                dihed_prms[ft][dihedtype] = [params]
            
            
        # Sections with parameters for the molecule
        #-------------------------------------------
        if section == 'moleculetype':
            molecule_exists=True
            molname = line.split()[0]
            if molname in molecules:
                print('\n%s'%('[ moleculetype ]'))
                print(line)
                # Update the shift to get the right index for atoms
                nshift = iat
                #print >> sys.stderr, 'nshift is %i when enter %s'%(nshift,molname)
            else:
                molecule_exists=False
        elif section == 'atoms' and molecule_exists:
            # Atom info was picked on the first loop
            # Get LJ parameters from database
            atoms[iat].setLJ(*atom_prms[atoms[iat].attype].split()[3:5])
            # Once LJ parameters are taken, we can change attype by attypeb (for OPLS)
            # WARNING: the link between atom_prms and the attype will be lost for OPLS at this point
            atoms[iat].attype = atom_prms[atoms[iat].attype].split()[-1]
            print(line) #+ '; mass:'+atom_prms[atoms[iat].attype].split()[0]
            # Update iat
            #print >> sys.stderr, 'Increasing iat to %i for after atom %s'%(iat+1,atoms[iat].attype)
            iat += 1
            
        elif section == 'bonds' and molecule_exists:
            comment=''
            data = line.split()
            i1=int(data[0])-1 + nshift
            i2=int(data[1])-1 + nshift
            ft=int(data[2])
            # if there is a dict for this ft. Otherwise, create to 
            # avoid errors
            if not ft in bond_prms:
                    bond_prms[ft] = dict()
            itemtype   = atoms[i1].attype+'-'+atoms[i2].attype
            itemtype_r = atoms[i2].attype+'-'+atoms[i1].attype
            if len(data) == 3 and ft != 5:
                # Get data from bondtype (consider also reverse order)
                itemtype_r = atoms[i2].attype+'-'+atoms[i1].attype
                if itemtype in bond_prms[ft]:
                    data.append(bond_prms[ft][itemtype])
                elif itemtype_r in bond_prms[ft]:
                    data.append(bond_prms[ft][itemtype_r])
                else:
                    data.append("0  0")
                    comment='; Not Found: '
                    #print '; '+itemtype
                    #print '; '+itemtype_r
                    #sys.exit()
            else:
                # Take the extra data as params
                params   = '   '.join(data[3:])
                del data[3:]
                #[ data.pop(-1) for i in range(len(data[3:]))]
                data.append(params)
            # Form bonds
            bonds.append(bond())    
            bonds[-1].setbond(*data)
            #bonds[-1].printbond()
            print(comment+bonds[-1].entryline()+'      ;'+itemtype)
            #----------------------
            #Print gaussian params
            #----------------------
            if itemtype not in gau_bonds and itemtype_r not in gau_bonds:
                gau_bonds.append(itemtype)
                if (bonds[-1].ft == 1):
                    r0_g = float(bonds[-1].prms.split()[0])*10.
                    kb_g = float(bonds[-1].prms.split()[1])/4.184/100.
                    print('HrmStr1 %s  %12.4f %12.4f'%(itemtype.replace('-','  '),kb_g,r0_g),file=f)
            #
            
        elif section == 'pairs' and molecule_exists:
            data = line.split()
            i1=int(data[0])-1 + nshift
            i2=int(data[1])-1 + nshift
            ft=int(data[2])
            # if there is a dict for this ft. Otherwise, create to 
            # avoid errors
            if not ft in pair_prms:
                    pair_prms[ft] = dict()
            itemtype   = atoms[i1].attype+'-'+atoms[i2].attype
            itemtype_r = atoms[i2].attype+'-'+atoms[i1].attype
            if len(data) == 3:
                # Get data from bondtype (consider also reverse order)
                itemtype_r = atoms[i2].attype+'-'+atoms[i1].attype
                if itemtype in pair_prms[ft]:
                    params = pair_prms[ft][itemtype]
                elif itemtype_r in pair_prms[ft]:
                    params = pair_prms[ft][itemtype_r]
                else:
                    # Generate pairs
                    # (V is sigma(cr=2,3) or C6(cr=1))
                    # (W is epsilon(cr=2,3) or C12(cr=1))
                    V1 = atoms[i1].V
                    W1 = atoms[i1].W
                    V2 = atoms[i2].V
                    W2 = atoms[i2].W
                    if combrule == 1:
                        C6  = np.sqrt(V1*V2) * fudgeLJ**(1./6.)
                        C12 = np.sqrt(W1*W2) * fudgeLJ**(1./12.)
                        params = "%12.4f %12.4f" % (C6,C12)
                    elif combrule == 2:
                        sgm = (V1+V2)/2.0
                        eps = np.sqrt(W1*W2) * fudgeLJ
                        params = "%12.4f %12.4f" % (sgm,eps)
                    elif combrule == 3:
                        sgm = np.sqrt(V1*V2) 
                        eps = np.sqrt(W1*W2) * fudgeLJ
                        params = "%12.4f %12.4f" % (sgm,eps)
                    else:
                        sys.exit('Unknown combrule: '+str(combrule))
            else:
                # Take the extra data as params
                params   = '   '.join(data[3:])
                del data[3:]
                
            # Change to ft=2
            if (ft==1):
                ft=2
                data[2]='2'
                if len(params.split()) == 2:
                    params = '%6.2f %8.3f %8.3f %s'%(fudgeQQ,atoms[i1].q,atoms[i2].q,params)
                
            # Add the parameters
            data.append(params)
                
            # Form pairs
            pairs.append(pair())    
            pairs[-1].setpair(*data)
            print(comment+pairs[-1].entryline()+'      ;'+itemtype)
            
        elif section == 'exclusions' and molecule_exists:
            print(line)
    
        elif section == 'angles' and molecule_exists:
            comment=''
            data = line.split()
            i1=int(data[0])-1 + nshift
            i2=int(data[1])-1 + nshift
            i3=int(data[2])-1 + nshift
            ft=int(data[3])
            itemtype   = atoms[i1].attype+'-'+atoms[i2].attype+'-'+atoms[i3].attype
            itemtype_r = atoms[i3].attype+'-'+atoms[i2].attype+'-'+atoms[i1].attype
            if len(data) == 4:
                # Get data from bondtype (consider also reverse order)
                itemtype_r = atoms[i3].attype+'-'+atoms[i2].attype+'-'+atoms[i1].attype
                if itemtype in angle_prms[ft]:
                    data.append(angle_prms[ft][itemtype])
                elif itemtype_r in angle_prms[ft]:
                    data.append(angle_prms[ft][itemtype_r])
                else:
                    data.append("0  0  0")
                    comment='; Not Found: '
            else:
                # Take the extra data as params
                params   = '   '.join(data[4:])
                del data[4:]
                data.append(params)
            # Form element
            angles.append(angle())    
            angles[-1].setangle(*data)
            print(comment+angles[-1].entryline()+'      ;'+itemtype)
            #----------------------
            #Print gaussian params
            #----------------------
            if itemtype not in gau_angles and itemtype_r not in gau_angles:
                gau_angles.append(itemtype)
                if (angles[-1].ft == 1):
                    a0_g = float(angles[-1].prms.split()[0])
                    ka_g = float(angles[-1].prms.split()[1])/4.184
                    print('HrmBnd1 %s  %12.4f %12.4f'%(itemtype.replace('-','  '),ka_g,a0_g),file=f)
            #
            
        elif section == 'dihedrals' and molecule_exists:
            comment=''
            data = line.split()
            i1=int(data[0])-1 + nshift
            i2=int(data[1])-1 + nshift
            i3=int(data[2])-1 + nshift
            i4=int(data[3])-1 + nshift
            ft=int(data[4])
            itemtype   = atoms[i1].attype+'-'+atoms[i2].attype+'-'+atoms[i3].attype+'-'+atoms[i4].attype
            if len(data) == 5:
                # Get data from bondtype. For dihedrals, we have a list of params (consider also reverse order)
                itemtype_r = atoms[i4].attype+'-'+atoms[i3].attype+'-'+atoms[i2].attype+'-'+atoms[i1].attype
                if itemtype in dihed_prms[ft]:
                    dihedral_param_list=dihed_prms[ft][itemtype]
                elif itemtype_r in dihed_prms[ft]:
                    dihedral_param_list=dihed_prms[ft][itemtype_r]
                # Try wild cards
                else:
                    # Search keys manually (only if have 'X' atoms)
                    Xitemtype = None
                    for dihed_key in dihed_prms[ft].keys():
                        dihed_bd = dihed_key.split('-')
                        if not 'X' in dihed_bd:
                            continue
                        dihed_fw = itemtype.split('-')
                        dihed_bw = itemtype_r.split('-')
                        # Try fw
                        found_params=True
                        for i in range(4):
                            if dihed_bd[i] != dihed_fw[i] and dihed_bd[i] != "X":
                                found_params=False
                        if found_params:
                            Xitemtype = dihed_key
                            break
                        # Try bw
                        found_params=True
                        for i in range(4):
                            if dihed_bd[i] != dihed_bw[i] and dihed_bd[i] != "X":
                                found_params=False
                        if found_params:
                            Xitemtype = dihed_key
                            break
                                
                    if Xitemtype in dihed_prms[ft]:
                        dihedral_param_list=dihed_prms[ft][Xitemtype]
                    else:
                        dihedral_param_list=["0  0  0"]
                        comment='; Not Found: '
                        #print itemtype
                        #print itemtype_r
                        #sys.exit()
            else:
                # Take the extra data as params
                params   = '   '.join(data[5:])
                del data[5:]
                dihedral_param_list=[params]
            # Form element(s)
            for i,params in enumerate(dihedral_param_list):
                # Remove comments (gromacs style)
                params = params.split(';')[0]
                # Remove comments (Joye style)
                params = params.split('#')[0]
                data.append(params)
                diheds.append(dihed())  
                diheds[-1].setdihed(*data)
                print(comment+diheds[-1].entryline()+'      ;'+itemtype+' # '+str(i+1))
                # Restore data for next cycle
                del data[-1]
                #----------------------
                #Print gaussian params
                #----------------------
                # Generate itemtype_and_params string with unique identifies for itemtype+params
                # so that two identical itemtypes with different params are printed
                plist_chr=[str(float(x)) for x in params.split()]
                itemtype_and_params='*'.join([itemtype]+plist_chr)
                itemtype_r_and_params='*'.join([itemtype_r]+plist_chr)
                if itemtype_and_params not in gau_diheds and itemtype_r_and_params not in gau_diheds:
                    gau_diheds.append(itemtype_and_params)
                    if (diheds[-1].ft == 1 or diheds[-1].ft == 9):
                        # Proper dihedral: Fourier type
                        p0_g = float(diheds[-1].prms.split()[0])
                        kd_g = float(diheds[-1].prms.split()[1])/4.184
                        n_g  = int(diheds[-1].prms.split()[2])
                        if (n_g > 4):
                            # This function type needs revision
                            print('! Warning: too large multiplicity for an AmbTrs. Maybe DreiTrs may fit:',file=f)
                            p0_g = p0_g + 180.
                            print('! DreiTrs %s  %12.4f %12.4f %2.1f %2.1f'%(itemtype.replace('-','  '),kd_g,p0_g,float(n_g),1.),file=f)
                        elif (n_g > 0):
                            if itemtype not in gau_diheds_amb.keys():
                                gau_diheds_amb[itemtype] = [[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0]]
                            gau_diheds_amb[itemtype][n_g-1] = [kd_g,p0_g]
                    elif (diheds[-1].ft == 2):
                        # Imnproper dihedral: harmonic => Transform to Fourier type (n=1)
                        # as harmonic types seems to be absent in Gaussian
                        # PO1 = phi_equil+180.0
                        # k1  = kharm
                        p0_g = float(diheds[-1].prms.split()[0])+180.
                        kd_g = float(diheds[-1].prms.split()[1])/4.184
                        n_g  = 1
                        print('ImpTrs %s  %12.4f %12.4f %2.1f'%(itemtype.replace('-','  '),kd_g,p0_g,float(n_g)),file=f)
                    elif (diheds[-1].ft == 3):
                        # RB type:
                        c0  = float(diheds[-1].prms.split()[0])
                        c1  = float(diheds[-1].prms.split()[1])
                        c2  = float(diheds[-1].prms.split()[2])
                        c3  = float(diheds[-1].prms.split()[3])
                        c4  = float(diheds[-1].prms.split()[4])
                        # Tranform to Fourier type first (see Gromacs manual)
                        k4_g = -c4/4.
                        k3_g = -c3/2.
                        k2_g = -c2 + 4.*k4_g
                        k1_g = -2*c1 + 3*k3_g
                        k4_g /= (2.*4.184)
                        k3_g /= (2.*4.184)
                        k2_g /= (2.*4.184)
                        k1_g /= (2.*4.184)
                        gau_diheds_amb[itemtype] = [[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0]]
                        gau_diheds_amb[itemtype][0] = [k1_g,0]
                        gau_diheds_amb[itemtype][1] = [k2_g,180]
                        gau_diheds_amb[itemtype][2] = [k3_g,0]
                        gau_diheds_amb[itemtype][3] = [k4_g,180]
                        
                #
            
        elif section == 'system':
            print(line)
            
        elif section == 'molecules':
            print(line)
            
    #print >> sys.stderr, ""
    #print >> sys.stderr, "Summary of gau_diheds"
    #for item in gau_diheds:
        #print >> sys.stderr, item
    #print >> sys.stderr, ""
    #print >> sys.stderr, "--End of summary"
    #print >> sys.stderr, ""

    # Print now Amber torisions for gaussian
    for dihed in gau_diheds_amb:
        p1 = int(gau_diheds_amb[dihed][0][1])
        p2 = int(gau_diheds_amb[dihed][1][1])
        p3 = int(gau_diheds_amb[dihed][2][1])
        p4 = int(gau_diheds_amb[dihed][3][1])
        v1 = gau_diheds_amb[dihed][0][0]
        v2 = gau_diheds_amb[dihed][1][0]
        v3 = gau_diheds_amb[dihed][2][0]
        v4 = gau_diheds_amb[dihed][3][0]
        print('AmbTrs %s %4i %4i %4i %4i %12.4f %12.4f %12.4f %12.4f %2.1f'%(dihed.replace('-','  '),\
                                                                             p1,p2,p3,p4,\
                                                                             v1,v2,v3,v4,1.),file=f)
        
    f.close()
