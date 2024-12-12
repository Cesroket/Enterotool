import pandas as pd
from rnanorm import CPM

def normalize_to_cpm(input_file, output_file):
    """
    Normaliza los datos de expresión génica en CPM utilizando rnanorm.

    :param input_file: Ruta al archivo TSV con los datos de expresión génica. Debe contener tres columnas: categoría, genotipo y cuentas.
    :param output_file: Ruta al archivo TSV donde se guardará el resultado normalizado.
    """
    # Leer el archivo de entrada
    print("Leyendo los datos de expresión...")
    exp_data = pd.read_csv(input_file, sep="\t", header=None, names=["Category", "Genotype", "Counts"])

    # Pivotar la tabla para tener genotipos como columnas y categorías como filas
    print("Reformateando los datos...")
    pivot_data = exp_data.pivot_table(index="Category", columns="Genotype", values="Counts", fill_value=0)

    # Normalizar a CPM
    print("Normalizando a CPM...")
    cpm = CPM().set_output(transform="pandas")
    normalized_data = cpm.fit_transform(pivot_data)

    # Guardar el resultado normalizado
    print("Guardando los datos normalizados en:", output_file)
    normalized_data.to_csv(output_file, sep="\t")
    print("Proceso completado.")

if __name__ == "__main__":
    import argparse

    # Definir argumentos para el script
    parser = argparse.ArgumentParser(description="Normalización de datos de expresión génica en CPM.")
    parser.add_argument("--input", required=True, help="Archivo TSV con los datos de expresión génica cruda.")
    parser.add_argument("--output", required=True, help="Archivo TSV donde se guardará el resultado normalizado.")

    # Parsear los argumentos
    args = parser.parse_args()

    # Ejecutar la normalización
    normalize_to_cpm(args.input, args.output)
