lidar-v2.1

Lidar Products for Elastic Scattering - SAVERNet Project

Este programa está destinado a obtener los productos de datos Lidar para los canales 
elásticos 1064 nm y 532 nm en sus dos componentes perpendiculares. Uno de los objetivos 
más importante de estos productos es caracterizar la presencia de aerosoles en la atmósfera, 
típicamente por medio del coeficiente de extinción o atenuación de aerosoles. No obstante, 
la atenuación debido a aerosoles no puede ser derivada directamente de la señal lidar de 
retrodispersión elástica, en tanto que es necesario asumir un valor de la razón lidar, 
i.e., la relación entre los coeficientes de extinción y retrodispersión de aerosoles, en 
los algoritmos de inversión. En consecuencia, el método utilizado aquí conduce a una importante 
incerteza en los productos finales. Por otro lado, la sencillez del método empleado lo hace 
apropiado para emplearse en una implementación operativa, tal como es nuestro propósito final, 
mientras que los resultados obtenidos pueden aún proveer información valiosa.

El métodos empleado está basado en los mismos principios empleados para obtener los 
productos operativos del National Institute for Environmental Studies, Japón. Una 
descripción de este método puede encontrarse en Shimizu et al. [2017]. Estos algoritmos 
han sido desarrollados para ser aplicados en la red de instrumentos lidares de AD-Net 
(Asian dust and aerosol lidar observation network) y han sido probados y validados 
principalmente para al caso  de polvo mineral asiático. Nosotros hemos implementado estos 
algoritmos con el fin de aplicarlos a la red de lidares ubicados en Argentina y Chile en 
el marco del proyecto SAVERNet (Sistema de Gestión de Riesgos Medioambientales Atmosféricos 
en Sudamérica). Si bien estos productos brindan información valiosa actualmente, aún se 
encuentran en fase experimental y se requiere al menos un año de pruebas y validaciones 
con datos externos independientes para adaptarlo a las condiciones de nuestra región.

Para ejecutarlo, usar:

./main.py Station

donde Station puede ser:
Aeroparque
Bariloche
Cordoba
Gallegos
Neuquen
Punta
Tucuman

Es importante definir todos los paths correctamente en el archivo de configuración:
parameters.cfg

Si se pretende realizar una ejecución automática para proveer los productos en tiempo 
casi real, se puede añadir la tarea:
make_all.sh
al Cron de Linux. Entonces puede ser necesario definir el path absoluto del programa 
en el encabezado del fichero:
main.py

Consultas y/o sugerencias:
lmingari@gmail.com
