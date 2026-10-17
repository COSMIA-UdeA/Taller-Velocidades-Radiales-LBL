# Taller: Extracción de Velocidades Radiales Utilizando el Método Line-by-Line

¡Bienvenido! En este repositorio encontrará material del taller realizado por el Semillero de Investigación en Astronomía COSMIA, liderado por estudiantes/egresados del Pregrado de Astronomía de la Universidad de Antioquia y egresados de la Universidad Pedagógica Nacional.

**¿Por qué hacer un taller de velocidades radiales?**

Los exoplanetas son uno de los temas de la astronomía más estudiados en la actualidad. Al día de hoy, conocemos más de 6000 exoplanetas confirmados y quedan muchos más por confirmar y ser descubiertos gracias a los nuevos surveys que saldrán en los próximos años. Descubrir exoplanetas alrededor de estrellas de tipo solar (K, F, G) se ha vuelto una tarea sistemática y más sencilla debido a los avances tecnológicos y herramientas matemáticas y computacionales disponibles en internet que analizan curvas de luz o datos de velocidad radial rápidamente entregando resultados de alta precisión. Por otro lado, es bien sabido que las estrellas más abundantes en nuestra galaxia son las enanas rojas o de tipo espectral M, que hospedan por lo menos un exoplaneta. Estas estrellas son magnéticamente muy activas, por lo que pueden provocar falsos positivos a la hora de detectar exoplanetas. Los métodos más utilizados para detectar exoplanetas a través de velocidades radiales (ej: template matching) no tienen en cuenta estos efectos. Por ende, COSMIA ha diseñado este taller con la intención de capacitar a los estudiantes de pregrado en teoría y algunas herramientas estándar sobre exoplanetas y velocidades radiales, y dotarlos de métodos a la vanguardia en la extracción de velocidades radiales para estrellas activas.

Para ello, utilizaremos el método Line-by-Line (ver Bouchy et al., 2001; Dumusque 2018; Cretignier et al., 2020; Artigau et al., 2022) que divide el espectro estelar en miles de líneas y calcula sus velocidades radiales para realizar posteriormente una suma ponderada obteniendo la velocidad radial final de la observación. Posteriormente, utilizaremos un modelo kepleriano con un proceso gaussiano para calcular los parámetros del sistema como la masa de los planetas. Estas técnicas las aplicaremos en algunas de las estrellas más conocidas y que más han causado debate sobre si son sistemas exoplanetarios, como Próxima Centauri o la estrella de Barnard.  

A continuación, encontrará un tutorial para la instalación de los paquetes [LBL](https://lbl.exoplanets.ca/) (Artigau et al., 2022) y [RadVel](https://radvel.readthedocs.io/en/latest/) (Fulton et al., 2018). Si bien estos paquetes tienen sus propios tutoriales de descarga y uso, creemos que no está de más una guía de instalación con lo necesario para el taller.

**Nota: Este tutorial ha sido diseñado para Linux. Asegúrate de tener instaladas y actualizadas las dependencias principales como Python, Pip, Git, etc.**

---
## Instalación: LBL
Para iniciar la instalación es necesario usar un environment (entorno virtual) como `conda` o `venv`. En este caso, utilizaremos este último.
> Nota del autor: ¿Qué es un virtual environment? Suponga que hace un tiempo instaló la versión 3.10 de Python en su sistema operativo. Sin embargo, para uno de sus proyectos se ve en la necesidad de utilizar un paquete que está optimizado para la versión 3.9. Para usted sería un problema instalar la versión 3.9 ya que, por una u otra razón, otro de sus proyectos corre únicamente con la versión 3.10. Las dos versiones no pueden coexistir al mismo tiempo en el mismo sistema. Por ende, crear un virtual environment le permite tener un "mini-sistema" dentro de su sistema donde únicamente se va a ejecutar la versión de Python de su preferencia sin afectar a la principal. Así, puede crear los environments que necesite para distintos proyectos. Es como dividir una casa en cuartos y que cada cuarto funcione de manera independiente. 

### Paso 1: Instalar las dependencias
Para usar `venv` es necesario instalar algunas librerías que Python necesita para compilarse a sí mismo.\
\
**Abre la terminal y corre:**
 
```bash
sudo apt update
sudo apt install -y \
  curl git \
  libbz2-dev libncurses-dev libffi-dev \
  libreadline-dev libssl-dev libsqlite3-dev \
  liblzma-dev zlib1g-dev tk-dev
```
_Este proceso solo lo necesitarás realizar una vez._

---
### Paso 2: Instalar `pyvenv`
**Instalar `pyenv` usando `curl`**
```bash
curl https://pyenv.run | bash
```
Bien, ahora se tuvo que haber creado una carpeta llamada `~/.pyenv`. Vamos a crear el archivo `.bashrc` que nos permitirá tener `pyenv` siempre disponible cuando abramos la terminal.
 
```bash
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
echo '[[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
echo 'eval "$(pyenv init - bash)"' >> ~/.bashrc
echo 'eval "$(pyenv virtualenv-init -)"' >> ~/.bashrc
```
> Nota del autor: ¿Para qué sirve cada línea? La primera línea apunta a la carpeta donde `~/.pyenv` fue instalado, la segunda línea permite que `pyenv` sea un comando reconocible, la tercera línea le permite a `pyenv` cambiar versiones y la cuarta línea activa una función extra de pyenv que no vamos a usar directamente, pero la incluimos para evitar mensajes de advertencia en la terminal.
 
**Ejecuta el archivo con `source`**
 
```bash
source ~/.bashrc
#pyenv --version   # Por si quieres verificar que todo salió bien
```
 
---

### Paso 3: Instalar Python 3.12

Acá vamos a instalar la versión 3.12.13

```bash
pyenv install 3.12.13
```
 
> Nota del autor: LBL necesita la versión 3.12. La razón para instalar la versión 3.12.13 es que el número de parche 13 es el más reciente para `pyenv` y podría evitar cualquier error. Sin embargo, el lector podría instalar la versión 3.12.0 y, en teoría, LBL correría sin problemas.
 
---
 
### Paso 4: Clonar el repositorio LBL

Dirígete a la carpeta desde donde vas a trabajar y clona el repositorio de LBL.
 
```bash
git clone https://github.com/njcuk9999/lbl.git
cd lbl
```
 
> Nota del autor: Artigau recomienda clonar el repositorio utilizando `git clone git@github.com:njcuk9999/lbl.git`. Sin embargo, encuentro más sencillo trabajar con el https.
 
---

### Paso 5: Crear el virtual enviroment
 
Es hora de crear el entorno. Para ello vamos a seguir la notación de Artigau "lbl-env", pero usted puede llamar a su entorno como desee.\
\
**Dentro de la carpeta `lbl`, configura la versión de Python y crea el entorno:**
 
```bash
pyenv local 3.12.13
python -m venv lbl-env
source lbl-env/bin/activate
```
 
La última línea activa el entorno. Debería aparecer `(lbl-env)` al inicio de tu línea de comando en la terminal.
 
> Nota del autor: `pyenv local` y `python -m venv` solo se necesitan correr una vez. En cambio, `source lbl-env/bin/activate` es el único comando que tú repites en cada sesión para activar tu entorno. Si quieres desactivarlo, solo escribe `deactivate`.
 
---
 
### Paso 6: Instalar LBL

Clonar el repositorio no nos instala LBL, para ello, con el entorno activo, instala LBL en modo desarrollador desde la carpeta clonada.
 
```bash
pip install -U -e .
```

---

## Instalación: RadVel 

No es necesario instalar RadVel desde el entorno creado, por lo que si quiere, puede desactivarlo y proceder con la instalación. Por fortuna, RadVel no necesita pasos complejos para descargarse, simplemente corra:

```bash
pip install radvel tables h5py
```
> Nota del autor: En el tutorial de RadVel ellos utilizan pytables en vez de tables. Esto causa un error porque pip no encuentra una dependencia con ese nombre. La razón para ello puede ser que utilizan conda o porque están utilizando el nombre con el cual el proyecto se conoce en la comunidad.

---

_Autor: @CamiloOspinal_\
_Editor:_\
_Versión:11.10.2026_\
_Un taller realizado por el equipo detrás de COSMIA_
