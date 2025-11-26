import json
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
from astroquery.vizier import Vizier
import numpy as np
import os
import configparser

class Utils:
    # Convertir coordenadas RA/DEC a objeto SkyCoord
    def obtener_coordenadas_estrella(ra_h, ra_m, ra_s, dec_h, dec_m, dec_s):
        ra = f"{ra_h}h{ra_m}m{ra_s}s"
        dec = f"{dec_h}d{dec_m}m{dec_s}s"
        coordenadas = SkyCoord(ra, dec, frame='icrs')
        return coordenadas

    # Calcular mínima elevación en el intervalo de tiempo
    def calcular_minima_elevacion(nombre_estrella, coordenadas_estrella, ubicacion, tiempos):
        altaz_frame = AltAz(obstime=tiempos, location=ubicacion)
        altaz_coordenadas = coordenadas_estrella.transform_to(altaz_frame)        
        elevaciones = altaz_coordenadas.alt.deg  # Elevaciones en grados
        elevacion_minima = np.min(elevaciones)        
        return elevacion_minima
     
    # Procesar JSON
    def procesar_archivo_json(ruta_archivo_json, output_json, ubicacion, tiempo_inicio, tiempo_fin, umbral_elevacion):
        # Cargar JSON
        with open(ruta_archivo_json, 'r') as archivo:
            datos = json.load(archivo)
        
        # Rango de tiempos
        tiempos = Time([tiempo_inicio, tiempo_fin])  
        
        # Nueva lista con solo las estrellas por encima el umbral de elevación mínima
        estrellas_filtradas = []
        total_estrellas = 0

        # Iterar sobre cada estrella en el JSON
        for estrella in datos['stars']:
            nombre = estrella['name']
            print(f"Estrella: {nombre}")

            # Obtener coordenadas RA/DEC
            coordenadas_estrella = Utils.obtener_coordenadas_estrella(
                estrella['ra_h'], estrella['ra_m'], estrella['ra_s'],
                estrella['dec_h'], estrella['dec_m'], estrella['dec_s']
            )

            # Calcular mínima elevación en el intervalo de tiempo
            elevacion_minima = Utils.calcular_minima_elevacion(nombre, coordenadas_estrella, ubicacion, tiempos)
            print(f"Elevación mínima de {nombre}: {round(elevacion_minima)} grados")

            # Filtrar estrellas con elevación mínima mayor al umbral
            if elevacion_minima > umbral_elevacion:
                #estrella['min_elevacion'] = elevacion_minima  # Agregar el dato de elevación minima(?)
                estrellas_filtradas.append(estrella)
                print(f"{elevacion_minima} > {umbral_elevacion} -> se guarda en {os.path.basename(output_json)} ")
            else:
                print(f"No se cumple {round(elevacion_minima)} > {umbral_elevacion} -> no se guarda en {os.path.basename(output_json)} ")
                
            total_estrellas = total_estrellas + 1
        
        # Se actualiza JSON con solo las estrellas que sobrepasan umbral de elevación mínima
        datos['stars'] = estrellas_filtradas

        # Se guarda JSON actualizado
        with open(output_json, 'w') as archivo_salida:
            json.dump(datos, archivo_salida, indent=4)

        print(f"Numero de estrella con elevación mínima por encima del limite {umbral_elevacion}: {len(estrellas_filtradas)},de un total de {total_estrellas}")
        print(f"Archivo JSON con estrellas filtradas guardado en: {output_json}")


# Ejemplo de uso
if __name__ == "__main__":
    
    dir_base = os.path.abspath( os.path.dirname( __file__ ))
    config = configparser.ConfigParser()

    # Leer config.ini
    config.read(os.path.abspath( os.path.dirname( __file__ )) + '/filtrar_json.ini')

    # Leer parámetros 
    input_json = os.path.abspath( os.path.dirname( __file__ )) + '/'+ config['RUTAS']['input_json']
    output_json = os.path.abspath( os.path.dirname( __file__ )) + '/'+ config['RUTAS']['output_json']
    lat = float(config['UBICACION']['latitud'])
    lon = float(config['UBICACION']['longitud'])
    altura = float(config['UBICACION']['altura'])
    ubicacion = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=altura*u.m)
    tiempo_inicio = config['TIEMPOS']['tiempo_inicio']
    tiempo_fin = config['TIEMPOS']['tiempo_fin']
    umbral_altura = int(config['SETTINGS']['umbral_altura'])

    # Valores leídos
    print(f"Archivo de entrada: {input_json}")
    print(f"Archivo de salida: {output_json}")
    print(f"Ubicación: {ubicacion}")
    print(f"Tiempo de inicio: {tiempo_inicio}")
    print(f"Tiempo de fin: {tiempo_fin}")
    print(f"Umbral: {umbral_altura}")      
    
    # Procesar archivo JSON
    Utils.procesar_archivo_json(input_json, output_json, ubicacion, tiempo_inicio, tiempo_fin, umbral_altura)

