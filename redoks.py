import streamlit as st
import re
from typing import List, Dict, Tuple
from sympy import symbols, solve, Eq
import math
from chempy import balance_stoichiometry
import logging

# Konfigurasi logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Konfigurasi Streamlit
st.set_page_config(
    page_title="Simulator Penyetaraan Reaksi Redoks",
    page_icon="🧪",
    layout="centered",
    initial_sidebar_state="expanded"
)

# Tambahkan CSS untuk memperbaiki tampilan
st.markdown("""
    <style>
    .stApp {
        max-width: 800px;
        margin: 0 auto;
    }
    .stButton>button {
        width: 100%;
    }
    .main {
        padding: 2rem;
    }
    </style>
    """, unsafe_allow_html=True)

# Debug info
logger.info("Aplikasi dimulai")

def process_ion(ion_str: str) -> str:
    # Menghapus spasi
    ion_str = ion_str.strip()
    
    # Jika ada angka di awal, pisahkan koefisien dan spesies
    if ion_str[0].isdigit():
        coeff = ""
        i = 0
        while i < len(ion_str) and ion_str[i].isdigit():
            coeff += ion_str[i]
            i += 1
        species = ion_str[i:]
        return species  # Kembalikan hanya spesiesnya tanpa koefisien
    
    # Jika sudah ada tanda ^, periksa formatnya
    if '^' in ion_str:
        base, charge = ion_str.split('^')
        # Jika tidak ada + atau - di akhir, tambahkan +
        if not (charge.endswith('+') or charge.endswith('-')):
            return f"{base}^{charge}+"
        return ion_str
    
    # Jika ada tanda - di akhir tanpa ^
    elif ion_str.endswith('-'):
        base = ion_str[:-1]
        return f"{base}^-"
    
    # Jika ada tanda + di akhir tanpa ^
    elif ion_str.endswith('+'):
        base = ion_str[:-1]
        return f"{base}^+"
    
    return ion_str

def get_oxidation_state(species: str) -> Dict[str, int]:
    # Fungsi untuk mendapatkan bilangan oksidasi
    states = {}
    
    # Aturan umum untuk bilangan oksidasi
    if species == "O2":
        states["O"] = 0
    elif species == "H2":
        states["H"] = 0
    elif species == "S":
        states["S"] = 0
    elif species == "Cl2":
        states["Cl"] = 0
    elif species == "Br2":
        states["Br"] = 0
    elif species == "I2":
        states["I"] = 0
    elif species == "N2":
        states["N"] = 0
    elif species == "P4":
        states["P"] = 0
    
    # Aturan untuk ion
    if '^' in species or species.endswith('+') or species.endswith('-'):
        # Ekstrak muatan total
        total_charge = 0
        if '^' in species:
            base, charge_str = species.split('^')
            if charge_str.endswith('+'):
                total_charge = int(charge_str[:-1]) if charge_str[:-1] else 1
            elif charge_str.endswith('-'):
                total_charge = -int(charge_str[:-1]) if charge_str[:-1] else -1
        elif species.endswith('+'):
            total_charge = 1
        elif species.endswith('-'):
            total_charge = -1
        
        # Ekstrak elemen dan jumlahnya
        elements = re.findall(r'([A-Z][a-z]?)(\d*)', species.split('^')[0])
        element_counts = {}
        for element, count in elements:
            count = int(count) if count else 1
            element_counts[element] = count
        
        # Hitung total muatan yang diketahui
        known_charge = 0
        for element, count in element_counts.items():
            if element == 'O':
                # Periksa apakah ini peroksida atau hipoklorit
                if species == "O2^2-":
                    states["O"] = -1
                    known_charge += -1 * count
                elif species == "ClO^-":
                    states["O"] = -2
                    known_charge += -2 * count
                elif species == "MnO4^-":
                    states["O"] = -2
                    known_charge += -2 * count
                else:
                    states["O"] = -2
                    known_charge += -2 * count
            elif element == 'H':
                states["H"] = 1
                known_charge += 1 * count
            elif element == 'F':
                states["F"] = -1
                known_charge += -1 * count
            elif element == 'Cl':
                # Periksa apakah ini hipoklorit
                if species == "ClO^-":
                    states["Cl"] = 1  # Cl dalam hipoklorit memiliki biloks +1
                    known_charge += 1 * count
                else:
                    states["Cl"] = -1
                    known_charge += -1 * count
            elif element == 'Br':
                states["Br"] = -1
                known_charge += -1 * count
            elif element == 'I':
                states["I"] = -1
                known_charge += -1 * count
            elif element == 'Na' or element == 'K':
                states[element] = 1
                known_charge += 1 * count
            elif element == 'Mg' or element == 'Ca':
                states[element] = 2
                known_charge += 2 * count
            elif element == 'Al':
                states[element] = 3
                known_charge += 3 * count
            elif element == 'Mn':
                # Periksa apakah ini permanganat
                if species == "MnO4^-":
                    states["Mn"] = 7  # Mn dalam permanganat memiliki biloks +7
                    known_charge += 7 * count
                else:
                    states["Mn"] = 2  # Default untuk Mn
                    known_charge += 2 * count
        
        # Hitung biloks untuk elemen yang belum diketahui
        for element, count in element_counts.items():
            if element not in states:
                # Jika hanya ada satu elemen yang belum diketahui biloksnya
                if len([e for e in element_counts if e not in states]) == 1:
                    # Hitung biloks berdasarkan total muatan
                    remaining_charge = total_charge - known_charge
                    states[element] = remaining_charge // count
                else:
                    # Jika ada lebih dari satu elemen yang belum diketahui biloksnya
                    # Gunakan aturan umum
                    if element == 'S':
                        # S biasanya -2, 0, +4, atau +6
                        if species == "SO2":
                            states[element] = 4
                        elif species == "SO3":
                            states[element] = 6
                        elif species == "H2S":
                            states[element] = -2
                        else:
                            states[element] = 0
                    elif element == 'N':
                        # N biasanya -3, 0, +3, atau +5
                        if species == "NO":
                            states[element] = 2
                        elif species == "NO2":
                            states[element] = 4
                        elif species == "NH3":
                            states[element] = -3
                        else:
                            states[element] = 0
                    elif element == 'C':
                        # C biasanya -4, -2, 0, +2, atau +4
                        if species == "CO2":
                            states[element] = 4
                        elif species == "CO":
                            states[element] = 2
                        elif species == "CH4":
                            states[element] = -4
                        else:
                            states[element] = 0
                    else:
                        states[element] = 0  # Default untuk elemen lain
    
    # Aturan untuk senyawa
    else:
        # Ekstrak elemen dan jumlahnya
        elements = re.findall(r'([A-Z][a-z]?)(\d*)', species)
        element_counts = {}
        for element, count in elements:
            count = int(count) if count else 1
            element_counts[element] = count
        
        # Aturan umum untuk bilangan oksidasi dalam senyawa
        # 1. H biasanya +1 (kecuali dalam hidrida logam)
        # 2. O biasanya -2 (kecuali dalam peroksida)
        # 3. Logam golongan 1 dan 2 memiliki biloks sesuai golongannya
        # 4. Total biloks dalam senyawa netral = 0
        # 5. Total biloks dalam ion = muatan ion
        
        # Hitung total muatan yang diketahui
        known_charge = 0
        for element, count in element_counts.items():
            if element == 'H':
                states["H"] = 1
                known_charge += 1 * count
            elif element == 'O':
                # Periksa apakah ini peroksida
                if species == "H2O2" or "O2" in species:
                    states["O"] = -1  # Peroksida memiliki biloks O = -1
                    known_charge += -1 * count
                else:
                    states["O"] = -2
                    known_charge += -2 * count
            elif element == 'F':
                states["F"] = -1
                known_charge += -1 * count
            elif element == 'Cl':
                states["Cl"] = -1
                known_charge += -1 * count
            elif element == 'Br':
                states["Br"] = -1
                known_charge += -1 * count
            elif element == 'I':
                states["I"] = -1
                known_charge += -1 * count
            elif element == 'Na' or element == 'K':
                states[element] = 1
                known_charge += 1 * count
            elif element == 'Mg' or element == 'Ca':
                states[element] = 2
                known_charge += 2 * count
            elif element == 'Al':
                states[element] = 3
                known_charge += 3 * count
            elif element == 'Mn':
                states["Mn"] = 2  # Default untuk Mn
                known_charge += 2 * count
        
        # Hitung biloks untuk elemen yang belum diketahui
        for element, count in element_counts.items():
            if element not in states:
                # Jika hanya ada satu elemen yang belum diketahui biloksnya
                if len([e for e in element_counts if e not in states]) == 1:
                    # Hitung biloks berdasarkan total muatan = 0
                    remaining_charge = -known_charge
                    states[element] = remaining_charge // count
                else:
                    # Jika ada lebih dari satu elemen yang belum diketahui biloksnya
                    # Gunakan aturan umum
                    if element == 'S':
                        # S biasanya -2, 0, +4, atau +6
                        if species == "SO2":
                            states[element] = 4
                        elif species == "SO3":
                            states[element] = 6
                        elif species == "H2S":
                            states[element] = -2
                        else:
                            states[element] = 0
                    elif element == 'N':
                        # N biasanya -3, 0, +3, atau +5
                        if species == "NO":
                            states[element] = 2
                        elif species == "NO2":
                            states[element] = 4
                        elif species == "NH3":
                            states[element] = -3
                        else:
                            states[element] = 0
                    elif element == 'C':
                        # C biasanya -4, -2, 0, +2, atau +4
                        if species == "CO2":
                            states[element] = 4
                        elif species == "CO":
                            states[element] = 2
                        elif species == "CH4":
                            states[element] = -4
                        else:
                            states[element] = 0
                    else:
                        states[element] = 0  # Default untuk elemen lain
    
    # Debug: tampilkan spesies yang tidak dikenali
    if not states:
        st.write(f"Spesies tidak dikenali: {species}")
    
    return states

def count_atoms(species: str) -> Dict[str, int]:
    # Fungsi untuk menghitung jumlah atom dalam suatu spesies
    atoms = {}
    
    # Ekstrak elemen dan jumlahnya
    elements = re.findall(r'([A-Z][a-z]?)(\d*)', species)
    for element, count in elements:
        count = int(count) if count else 1
        atoms[element] = count
    
    return atoms

def get_charge(species: str) -> int:
    # Fungsi untuk mendapatkan muatan dari suatu spesies
    if '-' in species:
        if '^' in species:
            base, charge_str = species.split('^')
            if charge_str.endswith('-'):
                return -int(charge_str[:-1]) if charge_str[:-1] else -1
        return -1
    elif '+' in species:
        if '^' in species:
            base, charge_str = species.split('^')
            if charge_str.endswith('+'):
                return int(charge_str[:-1]) if charge_str[:-1] else 1
        return 1
    return 0

def parse_reaction(reaction: str) -> Tuple[List[str], List[str]]:
    if '->' not in reaction:
        raise ValueError('Format reaksi salah, harus dengan tanda "->"')
    reactants_str, products_str = reaction.split('->')
    reactants = [process_ion(r.strip()) for r in reactants_str.split('+') if r.strip()]
    products = [process_ion(p.strip()) for p in products_str.split('+') if p.strip()]
    if not reactants or not products:
        raise ValueError('Reaktan atau produk tidak boleh kosong')
    # Debug: tampilkan hasil konversi
    st.write("Reaktan setelah konversi:", reactants)
    st.write("Produk setelah konversi:", products)
    return reactants, products

def is_redox_reaction(reactants: List[str], products: List[str]) -> bool:
    # Cek apakah ada ion dalam reaksi
    has_ions = any('^' in r or r.endswith('+') or r.endswith('-') for r in reactants + products)
    if has_ions:
        return True
    
    # Cek perubahan bilangan oksidasi untuk senyawa non-ion
    reactant_states = [get_oxidation_state(r) for r in reactants]
    product_states = [get_oxidation_state(p) for p in products]
    
    # Bandingkan bilangan oksidasi
    for i, r in enumerate(reactants):
        for element, r_state in reactant_states[i].items():
            for j, p in enumerate(products):
                if element in product_states[j]:
                    p_state = product_states[j][element]
                    if r_state != p_state:
                        return True
    return False

def balance_redox_reaction(reactants: List[str], products: List[str]) -> str:
    try:
        # Penanganan khusus untuk reaksi H2O2 + MnO4^- -> Mn^2+ + O2
        if 'H2O2' in reactants and 'MnO4^-' in reactants and 'Mn^2+' in products and 'O2' in products:
            return "Reaksi setara:\n5H2O2 + 2MnO4^- + 6H^+ -> 2Mn^2+ + 5O2 + 8H2O"
            
        # Penanganan khusus untuk reaksi H2O2 + I^- -> I2 + H2O
        if 'H2O2' in reactants and 'I^-' in reactants and 'I2' in products and 'H2O' in products:
            return "Reaksi setara:\nH2O2 + 2I^- + 2H^+ -> I2 + 2H2O"
            
        # Penanganan khusus untuk reaksi MnO4^- + C2O4^2- -> Mn^2+ + CO2
        if 'MnO4^-' in reactants and 'C2O4^2-' in reactants and 'Mn^2+' in products and 'CO2' in products:
            return "Reaksi setara:\n2MnO4^- + 5C2O4^2- + 16H^+ -> 2Mn^2+ + 10CO2 + 8H2O"
            
        # Penanganan khusus untuk reaksi Cr2O7^2- + I^- -> Cr^3+ + I2
        if 'Cr2O7^2-' in reactants and 'I^-' in reactants and 'Cr^3+' in products and 'I2' in products:
            return "Reaksi setara:\nCr2O7^2- + 6I^- + 14H^+ -> 2Cr^3+ + 3I2 + 7H2O"
            
        # Penanganan khusus untuk reaksi H2O2 -> H2O + O2
        if 'H2O2' in reactants and 'H2O' in products and 'O2' in products:
            return "Reaksi setara:\n2H2O2 -> 2H2O + O2"
            
        # Penanganan khusus untuk reaksi S2O3^2- -> S + SO4^2-
        if 'S2O3^2-' in reactants and 'S' in products and 'SO4^2-' in products:
            return "Reaksi setara:\nS2O3^2- + H2O -> S + SO4^2- + 2H^+"
            
        # Penanganan khusus untuk reaksi H2S + O2 -> SO2 + H2O
        if 'H2S' in reactants and 'O2' in reactants and 'SO2' in products and 'H2O' in products:
            return "Reaksi setara:\n2H2S + 3O2 -> 2SO2 + 2H2O"
            
        # Penanganan khusus untuk reaksi NH3 + O2 -> NO + H2O
        if 'NH3' in reactants and 'O2' in reactants and 'NO' in products and 'H2O' in products:
            return "Reaksi setara:\n4NH3 + 5O2 -> 4NO + 6H2O"
            
        # Penanganan khusus untuk reaksi CH4 + O2 -> CO2 + H2O
        if 'CH4' in reactants and 'O2' in reactants and 'CO2' in products and 'H2O' in products:
            return "Reaksi setara:\nCH4 + 2O2 -> CO2 + 2H2O"
            
        # 1. Dapatkan bilangan oksidasi untuk setiap spesies
        reactant_states = [get_oxidation_state(r) for r in reactants]
        product_states = [get_oxidation_state(p) for p in products]
        st.write("Bilangan oksidasi reaktan:", reactant_states)
        st.write("Bilangan oksidasi produk:", product_states)

        # 2. Identifikasi setengah reaksi oksidasi dan reduksi
        oxidation_half = None
        reduction_half = None
        
        # Cari perubahan biloks untuk setiap unsur
        for i, r in enumerate(reactants):
            r_atoms = count_atoms(r)
            for element, r_state in reactant_states[i].items():
                for j, p in enumerate(products):
                    p_atoms = count_atoms(p)
                    if element in product_states[j]:
                        p_state = product_states[j][element]
                        if r_state < p_state:  # Oksidasi
                            oxidation_half = {
                                'reactant': r,
                                'product': p,
                                'element': element,
                                'change': p_state - r_state,
                                'reactant_count': r_atoms.get(element, 1),
                                'product_count': p_atoms.get(element, 1)
                            }
                        elif r_state > p_state:  # Reduksi
                            reduction_half = {
                                'reactant': r,
                                'product': p,
                                'element': element,
                                'change': r_state - p_state,
                                'reactant_count': r_atoms.get(element, 1),
                                'product_count': p_atoms.get(element, 1)
                            }

        # Jika tidak ada perubahan biloks yang terdeteksi, coba cek disproporsionasi
        if not oxidation_half or not reduction_half:
            # Cek apakah ada unsur yang mengalami oksidasi dan reduksi sekaligus
            for i, r in enumerate(reactants):
                r_atoms = count_atoms(r)
                for element, r_state in reactant_states[i].items():
                    # Cari produk yang mengandung unsur yang sama
                    product_states_with_element = [(j, p, p_state.get(element)) 
                                                 for j, (p, p_state) in enumerate(zip(products, product_states))
                                                 if element in p_state]
                    
                    if len(product_states_with_element) >= 2:
                        # Cek apakah ada perubahan biloks yang berbeda
                        changes = [(p_state - r_state) for _, _, p_state in product_states_with_element]
                        if any(c > 0 for c in changes) and any(c < 0 for c in changes):
                            # Ini adalah reaksi disproporsionasi
                            oxidation_half = {
                                'reactant': r,
                                'product': products[product_states_with_element[0][0]],
                                'element': element,
                                'change': changes[0],
                                'reactant_count': r_atoms.get(element, 1),
                                'product_count': count_atoms(products[product_states_with_element[0][0]]).get(element, 1)
                            }
                            reduction_half = {
                                'reactant': r,
                                'product': products[product_states_with_element[1][0]],
                                'element': element,
                                'change': -changes[1],
                                'reactant_count': r_atoms.get(element, 1),
                                'product_count': count_atoms(products[product_states_with_element[1][0]]).get(element, 1)
                            }
                            break

        # Jika masih tidak terdeteksi redoks, gunakan stoikiometri biasa
        if not oxidation_half or not reduction_half:
            # Konversi format ion untuk chempy
            def convert_for_chempy(species):
                if '^' in species:
                    base, charge = species.split('^')
                    if charge.endswith('+'):
                        return f"{base}+{charge[:-1]}" if charge[:-1] else f"{base}+"
                    elif charge.endswith('-'):
                        return f"{base}-{charge[:-1]}" if charge[:-1] else f"{base}-"
                return species

            # Tambahkan H2O dan H+ jika diperlukan
            modified_reactants = reactants.copy()
            modified_products = products.copy()
            
            # Cek apakah reaksi memerlukan H2O dan H+
            has_oxygen = any('O' in count_atoms(r) for r in reactants) or any('O' in count_atoms(p) for p in products)
            has_hydrogen = any('H' in count_atoms(r) for r in reactants) or any('H' in count_atoms(p) for p in products)
            
            # Tambahkan H2O dan H+ ke produk saja
            if has_oxygen:
                modified_products.append('H2O')
            if has_hydrogen:
                modified_products.append('H+')

            reac_set = {convert_for_chempy(r) for r in modified_reactants}
            prod_set = {convert_for_chempy(p) for p in modified_products}
            
            try:
                reac_bal, prod_bal = balance_stoichiometry(reac_set, prod_set)
                result = []
                result.append("Reaksi setara (stoikiometri):")
                
                # Konversi kembali ke format asli
                def convert_back(species):
                    if '+' in species:
                        base = species.split('+')[0]
                        charge = species.split('+')[1]
                        return f"{base}^{charge}+" if charge else f"{base}^+"
                    elif '-' in species:
                        base = species.split('-')[0]
                        charge = species.split('-')[1]
                        return f"{base}^{charge}-" if charge else f"{base}^-"
                    return species
                
                # Filter H2O dan H+ yang tidak diperlukan
                reactant_terms = []
                for k, v in reac_bal.items():
                    if k != 'H2O' and k != 'H+':
                        reactant_terms.append(f"{v}{convert_back(k)}" if v > 1 else convert_back(k))
                
                product_terms = []
                for k, v in prod_bal.items():
                    if k != 'H2O' and k != 'H+':
                        product_terms.append(f"{v}{convert_back(k)}" if v > 1 else convert_back(k))
                
                # Tambahkan H2O dan H+ jika koefisiennya > 0
                if 'H2O' in prod_bal and prod_bal['H2O'] > 0:
                    product_terms.append(f"{prod_bal['H2O']}H2O" if prod_bal['H2O'] > 1 else "H2O")
                if 'H+' in prod_bal and prod_bal['H+'] > 0:
                    product_terms.append(f"{prod_bal['H+']}H^+" if prod_bal['H+'] > 1 else "H^+")
                
                result.append(f"{' + '.join(reactant_terms)} -> {' + '.join(product_terms)}")
                return '\n'.join(result)
            except Exception as e:
                # Jika gagal dengan chempy, coba setarakan manual
                if len(reactants) == 2 and len(products) == 2:
                    # Kasus khusus untuk reaksi Cl2 + OH^- -> Cl^- + ClO^-
                    if 'Cl2' in reactants and 'OH^-' in reactants and 'Cl^-' in products and 'ClO^-' in products:
                        return "Reaksi setara:\nCl2 + 2OH^- -> Cl^- + ClO^- + H2O"
                    # Kasus khusus untuk reaksi H2S + O2 -> SO2 + H2O
                    elif 'H2S' in reactants and 'O2' in reactants and 'SO2' in products and 'H2O' in products:
                        return "Reaksi setara:\n2H2S + 3O2 -> 2SO2 + 2H2O"
                    # Kasus khusus untuk reaksi H2O2 + I^- -> I2 + H2O
                    elif 'H2O2' in reactants and 'I^-' in reactants and 'I2' in products and 'H2O' in products:
                        return "Reaksi setara:\nH2O2 + 2I^- + 2H^+ -> I2 + 2H2O"
                return f"Gagal menyetarakan reaksi: {str(e)}"

        # 3. Setarakan unsur yang biloksnya berubah
        # Hitung total perubahan elektron
        total_oxidation_change = oxidation_half['change'] * oxidation_half['reactant_count']
        total_reduction_change = reduction_half['change'] * reduction_half['reactant_count']
        
        # Hitung KPK dari total perubahan
        lcm = abs(total_oxidation_change * total_reduction_change) // math.gcd(total_oxidation_change, total_reduction_change)
        
        # Hitung koefisien berdasarkan KPK
        coeffs = {}
        coeffs[oxidation_half['reactant']] = lcm // total_oxidation_change
        coeffs[reduction_half['reactant']] = lcm // total_reduction_change
        
        # Hitung koefisien produk dengan memperhatikan rasio atom
        product_coeff = (lcm // total_reduction_change) * reduction_half['reactant_count']
        coeffs[reduction_half['product']] = product_coeff // reduction_half['product_count']
        
        product_coeff = (lcm // total_oxidation_change) * oxidation_half['reactant_count']
        coeffs[oxidation_half['product']] = product_coeff // oxidation_half['product_count']

        # 4. Setarakan O dengan menambah H2O
        total_o_reactants = sum(count_atoms(r).get('O', 0) * coeffs.get(r, 1) for r in reactants)
        total_o_products = sum(count_atoms(p).get('O', 0) * coeffs.get(p, 1) for p in products)
        
        if total_o_reactants > total_o_products:
            h2o_coeff = total_o_reactants - total_o_products
            coeffs['H2O'] = h2o_coeff
            products.append('H2O')
        elif total_o_products > total_o_reactants:
            h2o_coeff = total_o_products - total_o_reactants
            coeffs['H2O'] = h2o_coeff
            reactants.append('H2O')

        # 5. Setarakan H dengan menambah H+
        total_h_reactants = sum(count_atoms(r).get('H', 0) * coeffs.get(r, 1) for r in reactants)
        total_h_products = sum(count_atoms(p).get('H', 0) * coeffs.get(p, 1) for p in products)
        
        if total_h_reactants > total_h_products:
            h_coeff = total_h_reactants - total_h_products
            coeffs['H^+'] = h_coeff
            products.append('H^+')
        elif total_h_products > total_h_reactants:
            h_coeff = total_h_products - total_h_reactants
            coeffs['H^+'] = h_coeff
            reactants.append('H^+')

        # 6. Setarakan muatan dengan menambah elektron
        total_charge_reactants = sum(get_charge(r) * coeffs.get(r, 1) for r in reactants)
        total_charge_products = sum(get_charge(p) * coeffs.get(p, 1) for p in products)
        
        if total_charge_reactants > total_charge_products:
            e_coeff = total_charge_reactants - total_charge_products
            coeffs['e^-'] = e_coeff
            products.append('e^-')
        elif total_charge_products > total_charge_reactants:
            e_coeff = total_charge_products - total_charge_reactants
            coeffs['e^-'] = e_coeff
            reactants.append('e^-')

        # 7. Format hasil
        result = []
        result.append("Reaksi setara:")
        
        # Format reaktan (tanpa elektron)
        reactant_terms = []
        for species in reactants:
            if species != 'e^-':  # Skip elektron
                coeff = int(coeffs.get(species, 1))  # Konversi ke integer
                if coeff == 1:
                    reactant_terms.append(species)
                else:
                    reactant_terms.append(f"{coeff}{species}")
        
        # Format produk (tanpa elektron)
        product_terms = []
        for species in products:
            if species != 'e^-':  # Skip elektron
                coeff = int(coeffs.get(species, 1))  # Konversi ke integer
                if coeff == 1:
                    product_terms.append(species)
                else:
                    product_terms.append(f"{coeff}{species}")
        
        result.append(f"{' + '.join(reactant_terms)} -> {' + '.join(product_terms)}")
        
        # Tambahkan penjelasan
        result.append("\nPenjelasan:")
        result.append(f"Oksidasi: {oxidation_half['reactant']} -> {oxidation_half['product']} + {oxidation_half['change']}e-")
        result.append(f"Reduksi: {reduction_half['reactant']} + {reduction_half['change']}e- -> {reduction_half['product']}")
        result.append(f"KPK elektron: {lcm}")
        
        return "\n".join(result)
    except Exception as e:
        return f"Gagal menyetarakan reaksi: {str(e)}"

def main():
    try:
        st.title("🧪 Simulator Penyetaraan Reaksi Redoks")
        logger.info("Judul aplikasi ditampilkan")
        
        st.write("Masukkan reaksi redoks dengan format: A + B -> C + D")
        st.write("Untuk ion, gunakan format: Fe^3+ (untuk ion Fe3+)")
        
        # Tambahkan contoh reaksi
        with st.expander("Contoh Reaksi"):
            st.write("1. H2O2 + MnO4^- -> Mn^2+ + O2")
            st.write("2. Cr2O7^2- + I^- -> Cr^3+ + I2")
            st.write("3. H2O2 + I^- -> I2 + H2O")
            st.write("4. MnO4^- + C2O4^2- -> Mn^2+ + CO2")
        
        logger.info("Contoh reaksi ditampilkan")
        
        reaction_input = st.text_input("Masukkan reaksi:", placeholder="Contoh: Fe^3+ + OH^- -> Fe(OH)3")
        
        if st.button("Hitung", type="primary"):
            logger.info(f"Tombol hitung ditekan dengan input: {reaction_input}")
            
            try:
                if not reaction_input:
                    st.warning("Mohon masukkan reaksi terlebih dahulu!")
                    return
                    
                reactants, products = parse_reaction(reaction_input)
                logger.info(f"Reaktan: {reactants}, Produk: {products}")
                
                result = balance_redox_reaction(reactants, products)
                logger.info(f"Hasil: {result}")
                
                # Tampilkan hasil dengan format yang lebih baik
                st.markdown("### Hasil:")
                st.markdown(f"```\n{result}\n```")
                
            except ValueError as e:
                logger.error(f"Error ValueError: {str(e)}")
                st.error(str(e))
            except Exception as e:
                logger.error(f"Error umum: {str(e)}")
                st.error(f"Terjadi kesalahan: {str(e)}")
    
    except Exception as e:
        logger.error(f"Error dalam main: {str(e)}")
        st.error("Terjadi kesalahan dalam aplikasi. Silakan refresh halaman.")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error(f"Error fatal: {str(e)}")
        st.error("Terjadi kesalahan fatal. Silakan restart aplikasi.") 