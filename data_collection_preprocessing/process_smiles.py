import os

from rdkit import Chem
# from standardiser import standardise
timeout = -1

_metal_nof = Chem.MolFromSmarts('[Li,Na,K,Rb,Cs,F,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]~[N,n,O,o,F]')
_metal_non = Chem.MolFromSmarts('[Al,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,Hf,Ta,W,Re,Os,Ir,Pt,Au]~[B,C,c,Si,P,As,Sb,S,Se,Te,Cl,Br,I,At]')
_metals = ['[Al]','[Sc]','[Ti]','[V]','[Cr]','[Mn]','[Fe]','[Co]','[Ni]','[Cu]','[Zn]','[Y]','[Zr]','[Nb]','[Mo]','[Tc]',
'[Ru]','[Rh]','[Pd]','[Pd++]','[Ag]','[Cd]','[Cd++]','[Hf]','[Ta]','[W]','[Re]','[Os]','[Ir]','[Pt]','[Au]','[Sn]','[Pb]','[Hg]','[Cd+2]','[Cr+3]','[Cr+6]',
'[Sn+3]','[Mg++]','[Sb+3]','[Al+3]','[Ba++]','[Fe+3]','[Mg+2]','[B+3]','[B]','[Pb++]','[Zr-2]','[Mn+2]','[Mn+3]','[Sb]','[Ti+4]','[Fe++]','[Pd++]',
'[Rb+]','[Mg]','[Ca++]','[As]','[Si]','[Ge]','[Te]','[Bi]','[Cu+]']

def disconnect(mol):
    """
    Adapated from molVS standardizer module. Now it returns the list of metals it has disconnected
    """
    metals = set([])
    for metal_atom in _metals:
        rwmol = Chem.RWMol(mol)
        smarts = Chem.MolFromSmarts(metal_atom)
        pairs = rwmol.GetSubstructMatches(smarts)
        for i, in reversed(pairs):
            metalSymbol = rwmol.GetAtomWithIdx(i).GetSmarts()
            metals.add(metalSymbol)
            rwmol.RemoveAtom(i)
        mol = rwmol.GetMol()

    for smarts in [_metal_nof, _metal_non]:
        pairs = mol.GetSubstructMatches(smarts)
        rwmol = Chem.RWMol(mol)
        orders = []
        for i, j in reversed(pairs):
            metalSymbol = mol.GetAtomWithIdx(i).GetSmarts()
            metals.add(metalSymbol)
            orders.append(int(mol.GetBondBetweenAtoms(i, j).GetBondTypeAsDouble()))
            rwmol.RemoveBond(i, j)
        # Adjust neighbouring charges accordingly
        mol = rwmol.GetMol()
        for n, (i, j) in enumerate(pairs):
            chg = orders[n]
            atom1 = mol.GetAtomWithIdx(i)
            atom1.SetFormalCharge(atom1.GetFormalCharge() + chg)
            atom2 = mol.GetAtomWithIdx(j)
            atom2.SetFormalCharge(atom2.GetFormalCharge() - chg)
    
    return mol, metals

