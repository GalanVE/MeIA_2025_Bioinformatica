with open('resultado600mol.txt', 'w') as arquivo:
    from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.Descriptors import HeavyAtomMolWt
from rdkit.Chem.Descriptors import FpDensityMorgan1
from rdkit.Chem.Descriptors import FpDensityMorgan2
from rdkit.Chem.Descriptors import FpDensityMorgan3
from rdkit.Chem.Descriptors import MaxPartialCharge
from rdkit.Chem.Descriptors import MaxAbsPartialCharge
from rdkit.Chem.Descriptors import MinPartialCharge
from rdkit.Chem.Descriptors import MinAbsPartialCharge
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.Descriptors import NumRadicalElectrons
from rdkit.Chem.Descriptors import NumValenceElectrons
from rdkit.Chem.Descriptors import BalabanJ
from rdkit.Chem.Descriptors import BertzCT
from rdkit.Chem.Descriptors import Ipc	
from rdkit.Chem.Descriptors import HallKierAlpha
from rdkit.Chem.Descriptors import Kappa1
from rdkit.Chem.Descriptors import Kappa3
from rdkit.Chem.Descriptors import Chi0
from rdkit.Chem.Descriptors import Chi1
from rdkit.Chem.Descriptors import Chi4n
from rdkit.Chem.Descriptors import Chi4v
from rdkit.Chem.Descriptors import MolLogP
from rdkit.Chem.Descriptors import MolMR
from rdkit.Chem.Descriptors import HeavyAtomCount	
from rdkit.Chem.Descriptors import NHOHCount
from rdkit.Chem.Descriptors import NOCount
from rdkit.Chem.Descriptors import NumHAcceptors
from rdkit.Chem.Descriptors import NumHDonors
from rdkit.Chem.Descriptors import NumHeteroatoms
from rdkit.Chem.Descriptors import NumRotatableBonds
from rdkit.Chem.Descriptors import NumAromaticRings
from rdkit.Chem.Descriptors import NumSaturatedRings
from rdkit.Chem.Descriptors import NumAliphaticRings
from rdkit.Chem.Descriptors import NumAromaticCarbocycles
from rdkit.Chem.Descriptors import NumSaturatedHeterocycles
from rdkit.Chem.Descriptors import NumSaturatedCarbocycles
from rdkit.Chem.Descriptors import NumAliphaticHeterocycles
from rdkit.Chem.Descriptors import NumAliphaticCarbocycles
from rdkit.Chem.Descriptors import RingCount
from rdkit.Chem.Descriptors import FractionCSP3
from rdkit.Chem.Descriptors import LabuteASA
from rdkit.Chem.Descriptors import PEOE_VSA1
from rdkit.Chem.Descriptors import PEOE_VSA14
from rdkit.Chem.Descriptors import SMR_VSA1
from rdkit.Chem.Descriptors import SMR_VSA10
from rdkit.Chem.Descriptors import VSA_EState1
from rdkit.Chem.Descriptors import VSA_EState10
from rdkit.Chem.Descriptors import Chi0n
from rdkit.Chem.Descriptors import Chi0v
from rdkit.Chem.Descriptors import EState_VSA1
from rdkit.Chem.Descriptors import EState_VSA11
from rdkit.Chem.rdMolDescriptors import CalcAUTOCORR2D
from rdkit.Chem.rdMolDescriptors import CalcNumSpiroAtoms
from rdkit.Chem.rdMolDescriptors import CalcCrippenDescriptors
from rdkit.Chem.rdMolDescriptors import CalcFractionCSP3
from rdkit.Chem.rdMolDescriptors import CalcHallKierAlpha
from rdkit.Chem.rdMolDescriptors import CalcLabuteASA
from rdkit.Chem.rdMolDescriptors import CalcNumAmideBonds
from rdkit.Chem.rdMolDescriptors import CalcNumUnspecifiedAtomStereoCenters
from rdkit.Chem.rdMolDescriptors import CalcNumAromaticHeterocycles
from rdkit.Chem.rdMolDescriptors import CalcNumBridgeheadAtoms
from rdkit.Chem.rdMolDescriptors import SlogP_VSA_
from rdkit.Chem.rdMolDescriptors import MQNs_


from csv import reader

read_smiles = open("/home/marcel/Documents/tabelaSmiles.csv",'r')

y = reader(read_smiles)

molecules = []

for y in read_smiles:
    molecules.append(y)




for mol in range(len(molecules)):
    
    print("Calculation SMILES: ", molecules[mol])
    m = Chem.MolFromSmiles(molecules[mol])
    print(Descriptors.TPSA(m),'TPSA')
    print(Descriptors.MolLogP(m),'LogP')
    AllChem.ComputeGasteigerCharges(m)
    print(float(m.GetAtomWithIdx(0).GetProp('_GasteigerCharge')),'Gasteiger')
    
    print(ExactMolWt(m),'mol weight')    
    print(HeavyAtomMolWt(m),'mol weight wht/ H')    
    print(FpDensityMorgan1(m),'Morgan 1')  
    print(FpDensityMorgan2(m),'Morgan 2')
    print(FpDensityMorgan3(m),'Morgan 3')
    print(MaxPartialCharge(m, force=False),'Max Carga Parcial')    
    print(MaxAbsPartialCharge(m, force=False),'Max ABS Carga Parcial')    
    print(MinPartialCharge(m, force=False),'Min Carga Parcial')    
    print(MinAbsPartialCharge(m, force=False),'Min ABS Carga Parcial')    
    print(MolWt(m),'mol weight total')    
    print(NumRadicalElectrons(m),'radical electron')    
    print(NumValenceElectrons(m),'valence electrons')
    
    
    print(BalabanJ(m),' BalabanJ')
    print(BertzCT(m),' BertzCT')
    print(Ipc(m),' Ipc')
    print(HallKierAlpha(m),'HallAlpha')
    print(Kappa1(m),'Kappa1')
    print(Kappa3(m),'Kappa3')
    print(Chi0(m),'Chi0')
    print(Chi1(m),'Chi1')
    print(Chi0n(m),'Chi0n')
    print(Chi4n(m),'Chi4n')
    print(Chi0v(m),'Chi0v' )
    print(Chi4v(m),'Chi4v')
    print(MolLogP(m),'MolLogP')
    print(MolMR(m),'MolMR')
    print(HeavyAtomCount(m),'HeavyAtomCount')
    print(NHOHCount(m),'NHOHCount')
    print(NOCount(m),'NOCount')
    print(NumHAcceptors(m),'NumHAcceptors')
    print(NumHDonors(m),'NumHDonors')
    print(NumHeteroatoms(m),'NumHeteroatoms')
    print(NumRotatableBonds(m),'NumRotatableBonds')
    print(NumAromaticRings(m),'NumAromaticRings')
    print(NumSaturatedRings(m),'NumSaturatedRings')
    print(NumAliphaticRings(m),'NumAliphaticRings')
    print(NumSaturatedHeterocycles(m),'NumSaturatedHeterocycles')
    print(NumSaturatedCarbocycles(m),'NumSaturatedCarbocycles')
    print(NumSaturatedRings(m),'NumSaturatedRings')
    print(NumAliphaticHeterocycles(m),'NumAliphaticHeterocycles')
    print(NumAliphaticCarbocycles(m),'NumAliphaticCarbocycles')
    print(NumAromaticCarbocycles(m),'NumAromaticCarbocycles')
    print(CalcNumAromaticHeterocycles(m),'CalcNumAromaticHeterocycles')
    print(CalcNumAmideBonds(m),'CalcNumAmideBonds')
    print(RingCount(m),'RingCount')
    print(FractionCSP3(m),'FractionCSP3')
    print(CalcNumBridgeheadAtoms(m),'CalcNumBridgeheadAtoms')
    print(LabuteASA(m),'LabuteASA')
    print(PEOE_VSA1(m),'PEOE_VSA1')
    print(PEOE_VSA14(m),'PEOE_VSA14')
    print(SMR_VSA1(m),'SMR_VSA1')
    print(SMR_VSA10(m),'SMR_VSA10')
    print(SlogP_VSA_(m),'SlogP_VSA_')
    print(EState_VSA1(m),'EState_VSA1')
    print(EState_VSA11(m),'EState_VSA11')
    print(VSA_EState1(m),'VSA_EState1')
    print(VSA_EState10(m),'VSA_EState10')
    print(MQNs_(m),'MQNs_')
    print(CalcAUTOCORR2D(m),'CalcAUTOCORR2D')
    print(CalcNumSpiroAtoms(m),'CalcNumSpiroAtoms')
    print(CalcCrippenDescriptors(m),'CalcCrippenDescriptors')
    print(CalcFractionCSP3(m),'CalcFractionCSP3')
    print(CalcHallKierAlpha(m),'CalcHallKierAlpha')
    print(CalcLabuteASA(m),'CalcLabuteASA')
    print(CalcNumUnspecifiedAtomStereoCenters(m),'CalcNumUnspecifiedAtomStereoCenters')
    print()
print('output600mol', file=arquivo)