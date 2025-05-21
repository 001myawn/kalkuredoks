import streamlit as st
import re
from typing import List, Dict, Tuple
from chempy import balance_stoichiometry

def parse_reaction(reaction: str) -> Tuple[List[str], List[str]]:
    if '->' not in reaction:
        raise ValueError('Format reaksi salah, harus dengan tanda "->"')
    reactants_str, products_str = reaction.split('->')
    reactants = [r.strip() for r in reactants_str.split('+') if r.strip()]
    products = [p.strip() for p in products_str.split('+') if p.strip()]
    if not reactants or not products:
        raise ValueError('Reaktan atau produk tidak boleh kosong')
    return reactants, products

def balance_redox_reaction(reactants: List[str], products: List[str]) -> str:
    try:
        reac_set = set(reactants)
        prod_set = set(products)
        reac_bal, prod_bal = balance_stoichiometry(reac_set, prod_set)
        output = []
        output.append("Reaksi setara:")
        left = " + ".join(f"{v if v > 1 else ''}{k}" for k, v in reac_bal.items())
        right = " + ".join(f"{v if v > 1 else ''}{k}" for k, v in prod_bal.items())
        output.append(f"{left} -> {right}")
        return "\n".join(output)
    except Exception as e:
        return f"Gagal menyetarakan reaksi otomatis: {e}"

def main():
    st.title("Simulator Penyetaraan Reaksi Redoks")
    st.write("Masukkan reaksi redoks dengan format: A + B -> C + D")
    reaction_input = st.text_input("Masukkan reaksi:", placeholder="Contoh: Fe + O2 -> Fe2O3")
    if st.button("Hitung"):
        try:
            reactants, products = parse_reaction(reaction_input)
            result = balance_redox_reaction(reactants, products)
            st.text(result)
        except ValueError as e:
            st.error(str(e))

default_example = "Fe + O2 -> Fe2O3"

if __name__ == "__main__":
    main() 