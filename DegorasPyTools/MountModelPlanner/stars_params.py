from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
import os
import json
import numpy as np
import math
import warnings
import traceback
import sys
warnings.simplefilter('ignore')

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

class Convert: 
    def name_to_fk5id(nombre_estrella):  # nombre -> fk5id
        
        try:
            result = sim.query_object(nombre_estrella) 
            
            if result is None:
                return f"No se encontró la estrella '{nombre_estrella}' en el catálogo SIMBAD."
            
            ids = result['IDS'][0]  # Se busca identificador FK5 en la lista de identificadores
            ids_list = ids.split('|')
            for identifier in ids_list:
                if 'FK5' in identifier:
                    return int(identifier.replace("FK5 ","")) 
                                
            return f"[Warning] la estrella '{nombre_estrella}' no tiene un id FK5 en SIMBAD."
                    
        except Exception as e:
            return f"Error en la consulta: {str(e)}"
            
    def fk5id_to_name(fk5id):     # fk5id --> nombre
        
        fk5id = "FK5 " + str(fk5id)
        
        try:
            # Hacemos una consulta por identificador
            result = sim.query_object(fk5id)
            
            if result is None:
                return f"No se encontró la estrella con el identificador FK5 '{fk5id}'."
            
            return result['MAIN_ID'][0]#.decode('utf-8')  # Se devuelve el nombre principal de la estrella
        
        except Exception as e:
            print(f"Error en la consulta: {str(e)}")

class Utils:

    def get_dict_star_param(ra_h, ra_m, ra_s, dec_h,  dec_m, dec_s, name, catalog_name, catalog_num, degoras_id, pm_ra, pm_dec, parallax, rad_vel):
        return {'ra_h': ra_h, 'ra_m': ra_m, 'ra_s': ra_s, 'dec_h': dec_h, 'dec_m': dec_m, 'dec_s': dec_s, 'name': name, 'catalog_name': catalog_name, 'catalog_num': catalog_num, 'degoras_id': degoras_id, 'pm_ra': pm_ra, 'pm_dec': pm_dec, 'parallax': parallax, 'rad_vel': rad_vel}

    def get_star_lists(fichero):

        lista_estrellas = []
        lista_estrellas_json = []
        
        with open(fichero) as fin:
            lineas = fin.readlines()
            for linea in lineas:
                item_fichero = linea.strip().split('#')[0].strip()
                if(item_fichero):
                    #print(f"-> {item_fichero}")
                    try:
                        fk5id = int(item_fichero)
                        nombre = Convert.fk5id_to_name(fk5id) 
                    except ValueError:
                        nombre = item_fichero
                        fk5id = Convert.name_to_fk5id(item_fichero) 
                    
                    if(fk5id != 0):
                        print(f"{nombre}, FK5id: {fk5id}")
                        if(not fk5id in lista_fk5id):
                            lista_estrellas.append(nombre)
                            lista_fk5id.append(fk5id)
                        else:
                            print(str(fk5id) + " ya ha aparecido, se omite")
                        
        return lista_estrellas, lista_fk5id

    def get_param_from_vizier_table(table_result):

        name = estrella
        catalog_num=table_result[0][0]
        ra_h = int(table_result[0][1].split(" ")[0])
        ra_m = int(table_result[0][1].split(" ")[1])
        ra_s = float(table_result[0][1].split(" ")[2])
        pm_ra = float(table_result[0][2])/100
        dec_h = int(table_result[0][3].split(" ")[0])
        dec_m = int(table_result[0][3].split(" ")[1])
        dec_s = float(table_result[0][3].split(" ")[2])
        pm_dec = float(table_result[0][4])/100
        parallax = float(table_result[0][5])
        if(math.isnan(parallax)):
            parallax = 0    
        rad_vel = float(table_result[0][6])
        degoras_id = catalog_num
        catalog_name = "FK5"
        
        return ra_h, ra_m, ra_s, dec_h,  dec_m, dec_s, name, catalog_name, catalog_num, degoras_id, pm_ra, pm_dec, parallax, rad_vel

"""
INICIO
"""

# Simbad
sim = Simbad()
sim.add_votable_fields('ids') 
# Vizier
viz = Vizier(columns=['FK5','RAJ2000', 'pmRA','DEJ2000','pmDE','plx', 'RV'],catalog='I/149A/catalog')

dir_base = os.path.abspath( os.path.dirname( __file__ ))  # Carpeta donde está este script
lista_estrellas = []
lista_fk5id = []
lista_estrellas_json = []
num_estrellas = 0

n_arg = len(sys.argv) # Argumentos linea de comandos

print("***\n1.- Identificar estrellas")
if(n_arg == 1):
    lista_estrellas, lista_fk5id = Utils.get_star_lists(dir_base + "/input.txt")
else:

    try:
        fk5id = int(sys.argv[1])
        nombre = Convert.fk5id_to_name(fk5id) 
    except ValueError:
        nombre = sys.argv[1]
        fk5id = Convert.name_to_fk5id(sys.argv[1]) 

    lista_estrellas, lista_fk5id = [nombre],[fk5id]

print("***\n2.- Obtener parametros de las estrellas")

for estrella in lista_estrellas:  
    
    try:
        #print(f"****\n{estrella}:")
        table_result = viz.query_object(estrella)['I/149A/catalog']
        #print(table_result)
        ra_h, ra_m, ra_s, dec_h,  dec_m, dec_s, name, catalog_name, catalog_num, degoras_id, pm_ra, pm_dec, parallax, rad_vel = Utils.get_param_from_vizier_table(table_result)
        lista_estrellas_json.append(Utils.get_dict_star_param(ra_h, ra_m, ra_s, dec_h,  dec_m, dec_s, name, catalog_name, catalog_num, degoras_id, pm_ra, pm_dec, parallax, rad_vel)  )
        num_estrellas = num_estrellas + 1
    except:
        #print(traceback.format_exc())
        print(f"Error en la consulta ({estrella}), se omite")

print("***\n3.- Crear fichero json")

jsondata = { "stars": lista_estrellas_json }
json_object = json.dumps(jsondata, indent = 4,cls=NpEncoder)  

nombre_fichero = "stars_params.json"
ruta_fichero = dir_base + "/" + nombre_fichero
with open(ruta_fichero, 'w', encoding='utf-8') as f:
    json.dump(jsondata, f, ensure_ascii=False, indent=4,cls=NpEncoder)
    print(f"Creado fichero: {nombre_fichero}")

print(f"****\nNúmero de estrellas procesadas correctamente: {num_estrellas}")


