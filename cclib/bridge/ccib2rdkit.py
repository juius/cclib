"""Bridge between cclib data and rdkit (https://www.rdkit.org/)."""

from packaging import version
from cclib.parser.data import ccData
from cclib.parser.utils import find_package

_found_rdkit = find_package("rdkit")
if _found_rdkit:
    import rdkit
    from rdkit import Chem
    _supports_xyz2mol = version.parse(rdkit.__version__) > version.parse('2022.09')
    if _supports_xyz2mol:
        from rdkit.Chem import rdDetermineBonds
    else:
        raise Warning(f'You use rdkit version {rdkit.__version__} which does not support the automatic \
            determination of bond orders and connectivity. Use version > 2022.9 for this feature.')
    
def _check_rdkit(_found_rdkit):
    if not _found_rdkit:
        raise ImportError("You must install `rdkit` to use this function")
    
def makecclib(mol):
    """Create cclib attributes and return a ccData from an rdkit molecule.

    Beyond the numbers, masses and coordinates, we could also set the total charge
    and multiplicity, but often these are calculated from atomic formal charges
    so it is better to assume that would not be correct.
    """
    _check_rdkit(_found_rdkit)
    attributes = {
        'atomcoords':   [],
        'atommasses':   [],
        'atomnos':      [],
        'natom':        mol.GetNumAtoms(),
    }
       
    for atom in mol.GetAtoms():
        attributes['atommasses'].append(atom.GetMass())
        attributes['atomnos'].append(atom.GetAtomicNum())

    for conf in mol.GetConformers():
        attributes['atomcoords'].append(conf.GetPositions())
        
    return ccData(attributes)

def makerdkit(data, firstgeom=False, lastgeom=True, allgeom=False):
    """Create an rdkit molecule."""
    from cclib.io import ccwrite
    _check_rdkit(_found_rdkit)
    xyz_string = ccwrite(data, 'xyz', firstgeom=firstgeom, lastgeom=lastgeom, allgeom=allgeom)
    rdkit_mol = Chem.MolFromXYZBlock(xyz_string)
    if _supports_xyz2mol:
        try:
            rdDetermineBonds.DetermineBonds(rdkit_mol, charge=data.charge)
        except ValueError as e:
            rdDetermineBonds.DetermineConnectivity(rdkit_mol, charge=data.charge)
            # this only one unsupported element
            print(f'Bond order can not be perceived since {str(e).split()[-1]} is not supported.') 
        
    return rdkit_mol

def copy_bond_order(mol, template_mol): 
    # assert all same atomtype between molecules
    assert [a.GetSymbol() for a in mol.GetAtoms()] == [a.GetSymbol() for a in template_mol.GetAtoms()],\
        'Missmatch in atoms between mol and template_mol.'
    # assert that template_mol contains bonds
    assert len(template_mol.GetBonds()) > 0, 'No bonds in template_mol'
    
    emol = Chem.EditableMol(mol)
    # Remove all bonds in mol
    _ = [emol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]
    # copy all bonds from template_mol
    for bond in template_mol.GetBonds():
        emol.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), order=bond.GetBondType())

    return emol.GetMol()