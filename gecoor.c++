/*
 * ****gecoor****
 * (Generador de coordenadas atómicas de nanotubos armchair con
 * longitudes y ángulos de enlace específicos.)
 * 
 * 
 * Autores:
 * 	Gustavo Dominguez Rodriguez[1]
 * 	Gabriel Ivan Canto Santana [1]
 * 	Jorge Alejandro Tapia Gonzalez [2]
 * 	Cesar Alberto Cab Cauich[2]
 * 
 * 	[1]	Centro de Investigacion en Corrosion, Universidad Autonoma de
 * 		Campeche.
 * 	[2]	Facultad de Ingenieria, Universidad Autonoma de Yucatan
 * 
 * 
 * 
 * ****Ecuaciones para la determinación de las posiciones atómicas
 * modificadas.****
 * 
 * El paquete considera cuatro parámetros de entrada, n, r_a, r_b y
 * θ_bb, los cuales son el índice de quiralidad en un nanotubo tipo
 * armchair de la forma (n,n), la distancia de los enlaces C-C
 * horizontales, las distancias de los enlaces C-C oblicuos y el ángulo
 * entre dos enlaces  C-C oblicuos. Con ellos se desea construir un
 * nanotubo completo.
 * 
 * Una manera sencilla es teniendo el índice de la quiralidad (n), la
 * altura de la celda unitaria (z), la distancia del centro del nanotubo
 * al punto medio de cada enlace horizontal (A), y la longitud de dichos
 * enlaces (r_a). De manera inmediata se cuenta con n y r_a. Por lo cual
 * deben hallar z y A.
 * 
 * Para ello se toma en cuenta que los centros de los enlaces
 * horizontales forman un polígono regular con un número total de lados
 * 2n. Ya que por cada n en (n,n) hay una cresta, que incluye un enlace
 * horizontal en posición vertical baja y otro horizontal en posición
 * vertical alta. La longitud de los lados L del polígono se encuentra
 * definido como definidos como:
 * 
 * 		L = 2A sin⁡(π/2n)
 * 
 * El polígono regular puede interpretarse como una secuencia de
 * triángulos isósceles con base igual a L, cuyo ángulo en la arista
 * compartida por todos se encuentra definido como
 * 
 * 		α = 2π/2n = π/n
 * 
 * Con dicho ángulo puede obtenerse el ángulo entre dos caras del
 * polígono regular (β). Sabiendo que los ángulos de un triángulo suman
 * 180 °, y que ambos ángulos opuestos a α son un medio
 * 
 *		π = α+2(β/2) = α+β
 * 		β = π-α = π-π/n,
 * 
 * Por otro lado, el ángulo formado entre un lado del polígono y el
 * enlace horizontal (γ) se puede obtener sabiendo que hay dos y son
 * complementarios con β.
 * 
 * 		π = β+2γ
 * 		γ = (π-β)/2 = (π-π+π/n)/2 = π/2n
 * 
 * Con ellos, ya se puede obtener la proyección (p_a) de los enlaces r_a
 * sobre L. Sabiendo que debido al punto de anclaje de la arista con el
 * enlace C-C está en el centro del mismo, r_a con tribuye con un medio
 * de su valor.
 * 
 * 		p_a = r_a/2 cos⁡(γ) = r_a/2 cos⁡(π/2n)
 * 
 * Sabiendo que la proyección de r_b en el plano transversal (p_b) es
 * paralela a los lados del polígono. Se puede obtener su valor
 * substrayendo dos veces la longitud de p_a a la longitud de las caras
 * (L)
 * 
 * 		p_b = L-2p_a = 2A sin⁡(π/2n)-2[r_a/2 cos⁡(π/2n)]
 * 			= 2A sin⁡(π/2n) - r_a cos⁡(π/2n)							(0)
 * 
 * Por otro lado, se sabe que r_b es la hipotenusa de un rectángulo en
 * conjunto con los catetos p_b y z/2,
 * 
 * 		r^2_b = p^2_b+z^2/4											(1)
 * 
 * Cuyo angulo entre r_b y p_b, equivalente a un medio de θ_(bb), podría 
 * obtenerse como:
 * 
 * 		θ_bb/2 = ArcTan[(z/2)/p_b] = ArcTan[z/(2p_b)]				(2)
 * 		Tan(θ_bb/2) = z/(2p_b)
 * 
 * Despejando z
 * 
 * 		z = 2p_b Tan(θ_bb/2)										(3)
 * 
 * Sustituyendo (3) en (1) se puede obtener el valor de la proyección
 * p_b en función de θ_bb
 * 
 * 		r^2_b = p^2_b+[2p_b Tan(θ_bb/2)]^2/4 = p^2_b+p^2_b Tan^2(θ_bb/2)
 * 			  = p^2_b[1+Tan^2(θ_bb/2)] = p^2_b/Cos^2(θ_bb/2)
 * 		p_b=r_b Cos(θ_bb/2)											(4)
 * 
 * Finalmente, para obtener z se sustituye, (4) en (3),
 * 
 * 		z = 2r_b Cos(θ_bb/2) Tan(θ_bb/2) = 2r_b Sin(θ_bb/2)
 * Mientras que para obtener A, se sustituye (4) en (0),
 * 
 * 		r_b Cos(θ_bb/2) = 2A sin⁡(π/2n)-r_a cos⁡(π/2n)
 *		A = [r_a cos⁡(π/2n)+r_b Cos(θ_bb/2)]/[2 sin⁡(π/2n)]
 * 
 * Por lo tanto, ya se conocen los cuatro parámetros requeridos para
 * construir un nanotubo que satisfaga las condiciones geométricas
 * requeridas.
 * 
 * 
 * 
 * 
 * ****Instalación y Uso.****
 * 
 * Para compilar el paquete basta con realizar la siguiente instrucción
 * en el directorio donde se encuentre el código fuente:
 * 
 * 		g++ -o gecoor -Wall gecoor.cpp -static-libgcc -static-libstdc++
 * 
 * lo que generará el archivo ejecutable gecoor. En algunas ocasiones es
 * necesario cambiarle los privilegios a ejecutable mediante la
 * instrucción
 * 
 * 		chmod u+x gecoor
 * 
 * Ya con ello, se podrá ejecutar por cualquier usuario. Siempre y
 * cuando se especifique la ruta completa.
 * 
 * Si no se desea especificar la ruta cada vez, puede agregarse a la
 * variable de entorno PATH, añadiendo el siguiente comando en el
 * archivo ".bashrc" del usuario
 * 
 * 		export PATH=$PATH:[ruta de instalación]/
 * 
 * Donde [ruta de instalación] es la ruta de donde se realizó la}
 * instalación, sin poner corchetes y empleando el carácter "/" como
 * separador de directorios.
 * 
 * Por otro lado, la ejecución del paquete gecoor es la siguiente:
 * 
 * 		gecoor [n] [U] [r_a] [r_b] [th_bb] 
 * 
 * donde n es el índice de quiralidad del nanotube armchair (n,n) , U es
 * el ancho de la celda unitaria en los componentes transversales en
 * Angstrongs, r_a es la longitud de los enlaces horizontales en
 * Angstrongs, r_b es la longitud de los enlaces oblicuos en Angstrongs,
 * y th_bb es el ángulo entre enlaces oblicuos en radianes. Cuando se
 * ejecuta genera imprime en pantalla el contenido de un archivo de
 * entrada para VASP tipo CONTCAR/POSCAR. Si el archivo se desea
 * almacenar en un archivo POSCAR vasta con redirigir la salida estándar
 * de la siguiente manera:
 * 
 * 		gecoor [n] [U] [r_a] [r_b] [th_bb] >POSCAR

*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;

float distancia(float *v1, float *v2){
	float x = v1[0]-v2[0];
	x *= x;
	float y = v1[1]-v2[1];
	y *= y;
	float z = v1[2]-v2[2];
	z *= z;

	return sqrt(x+y+z);
}

float distanciacentral(float *v1, float *v2){
	float x = 0.5*(v1[0]+v2[0]);
	x *= x;
	float y = 0.5*(v1[1]+v2[1]);
	y *= y;
	float z = 0.5*(v1[2]+v2[2]);
	z *= z;
	return sqrt(x+y+z);
}

float angulo(float *v1, float *v2, float *v3){
	float ax = v1[0] - v2[0], ay = v1[1] - v2[1], az = v1[2] - v2[2];
	float bx = v3[0] - v2[0], by = v3[1] - v2[1], bz = v3[2] - v2[2];

	return acos((ax*bx+ay*by+az*bz)/sqrt((ax*ax + ay*ay + az*az)*(bx*bx + by*by + bz*bz)));
}

float angulosimetricoz(float *v1, float *v2){
	float ax = v1[0] - v2[0], ay = v1[1] - v2[1], az = v1[2] - v2[2];
	float axx = ax*ax, ayy = ay*ay, azz = az*az;

	return acos((axx+ayy-azz)/(axx + ayy + azz));
}

int main(int argc, const char **args){
	if(argc != 6){
		cout <<"\n\nSintaxis:\n\gecoor <quiralidad> <ancho de celda> <r_a> <r_b> <th_bb>\n\n";

		return 1;
	}
	int n = atoi(args[1]);
	float U = atof(args[2]);
	float ra = atof(args[3]);
	float rb = atof(args[4]);
	float th = atof(args[5]);

	int N = 2*n;
	int num = 4*n;
	float th_s2 = 0.5*th;
	float alph = 3.141599265359/N;

	float z = 2.*rb*sin(th_s2);
	float r = (ra*cos(alph)+rb*cos(th_s2))/(2.*sin(alph));
	
	cout <<setprecision(17)  <<"CNT (" <<n <<"," <<n <<"): Ancho de celda (" <<U <<"), r_a(" <<ra <<"), r_b(" <<rb <<"), th_b-b(" <<th <<")\n";

	float bases[]={U, 0., 0.
		, 0., U, 0.
		, 0., 0., z};

	cout <<"   1.000000\n     "
		<<bases[0] <<"\t" <<bases[1] <<"\t" <<bases[2] <<"\n     "
		<<bases[3] <<"\t" <<bases[4] <<"\t" <<bases[5] <<"\n     "
		<<bases[6] <<"\t" <<bases[7] <<"\t" <<bases[8] <<"\n   C\n     " <<num <<"\nDirect\n";

	float **atomos;
	
	float rnorm = r/U;
	float ras2norm = ra/(2.*U);
	atomos = new float*[num];
	for(int i = 0; i < N; i++){
		atomos[2*i] = new float[3];
		atomos[2*i+1] = new float[3];
		
		float cai=cos(2.*alph*i);
		float sai=sin(2.*alph*i);

		float rac=ras2norm*cai;
		float ras=-ras2norm*sai;
		atomos[2*i][0] = rnorm*cai;
		atomos[2*i][1] = rnorm*sai;
		atomos[2*i][2] = 0.5*(i % 2);
		atomos[2*i+1][0] = atomos[2*i][0] + ras;
		atomos[2*i+1][1] = atomos[2*i][1] + rac;
		atomos[2*i+1][2] = atomos[2*i][2];
		atomos[2*i][0] -= ras;
		atomos[2*i][1] -= rac;
	}

	for(int i = 0; i < num; i++){
		cout <<"    " <<atomos[i][0] <<"\t" <<atomos[i][1] <<"\t" <<atomos[i][2] <<"\n";
		delete[] atomos[i];
	}
	delete[] atomos;
	
	cout <<"\n";

	for(int i = 0; i < num; i++){
		cout <<"    0.000000\t0.000000\t0.000000\n";
	}

	return 0;
}
