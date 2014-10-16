/*
CSCI 480
Assignment 3 Raytracer

Name: <Your name here>
*/
#include <iostream>
#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <pic.h>
#include <string.h>
#include <math.h>
#include <vector>
using namespace std;

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10
#define PI 3.14159265

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Point 
{
    double x;
    double y;
    double z;
    Point()
    {
    	x=0;y=0;z=0;
    }
    Point (double _x, double _y, double _z)
    {
    	x=_x;y=_y;z=_z;
    }
};

struct Ray 
{
	Point p;
	Point d;
};

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

struct Near
{
	double z,lightToObject,lightToSurface;
	Point color;
	//storage function
	Point bary,intersection,normal;
	Triangle triangle;
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];
bool antialiasing=false,softshadow=false;

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

Point subtract (Point p1, Point p2)
{
	Point p3;
	p3.x = p1.x - p2.x;
	p3.y = p1.y - p2.y;
	p3.z = p1.z - p2.z;
	return p3;
}

Point Cross(Point p1, Point p2) //Cross product
{
    Point product;
    product.x = p1.y*p2.z-p1.z*p2.y;
    product.y = p1.z*p2.x-p1.x*p2.z;
    product.z = p1.x*p2.y-p1.y*p2.x;
    return product;
}

float Dot(Point p1, Point p2) //Dot product
{
    return p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
}

Point Normalize (Point &p) //Normalize a vector
{
    float length = sqrt(Dot(p,p));
    p.x /= length; p.y /= length; p.z /= length;
    return p;
}

void Clamp (double &d)
{
	if (d < 0)
	{
		d = 0;
	}
	else if (d > 1)
	{
		d = 1;
	}
}

//MODIFY THIS FUNCTION
void draw_scene()
{
  unsigned int x,y;
  
  double a = (double)WIDTH/HEIGHT;
  Point upLeft,upRight,downLeft,downRight;
  upLeft.x = -a*tan(fov*PI/360); upLeft.y = tan(fov*PI/360); upLeft.z = -1;
  upRight.x = a*tan(fov*PI/360); upRight.y = tan(fov*PI/360); upRight.z = -1;
  downLeft.x = -a*tan(fov*PI/360); downLeft.y = -tan(fov*PI/360); downLeft.z = -1;
  downRight.x = a*tan(fov*PI/360); downRight.y = -tan(fov*PI/360); downRight.z = -1;
  
  Triangle triangle,tri;
  Sphere sphere,_sphere;
  Ray ray, shadowRay;
  Point pixel,colorTriangle,colorSphere,barycentric,picked,bary;
  double w = upRight.x - upLeft.x, h = upLeft.y - downLeft.y,qS,fS,uS,vS;
  double pixelWidth = w/WIDTH, pixelHeight = h/HEIGHT,Shininess,spec,_d,_areaABC,_areaPBC,_areaPCA,t5;
  double pixelCenterX = pixelWidth/2, pixelCenterY = pixelHeight/2,nx,ny,nz,redDiffuse,redSpecular,greenDiffuse,greenSpecular,blueDiffuse,blueSpecular;
  double pixelX, pixelY,b,c,t0,t1,t2,t3,t4,diffuse,specular,scalar,red,green,blue,_b,_c,_t0,_t1,d,i1,i2,i3,zTriangle,zSphere,areaABC,areaPBC,areaPCA,min;
  bool shadow,triangleB,sphereB;
  vector<Near> store;
  vector<Near> closest;
  vector<Point> softShadow;
  Near near,nearest;
  
  for(x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    
    pixelX = upLeft.x+(pixelWidth*x)+pixelCenterX;
    for(y=0;y < HEIGHT;y++)
    {
      pixelY = upLeft.y-(pixelHeight*y)-pixelCenterY;
      Point pixel (pixelX,pixelY,-1);
      ray.d = pixel;
      Normalize(ray.d);
      triangleB = false;
      sphereB = false;
      zTriangle = 0;
      colorTriangle = Point(0,0,0);
      int shadows=0;
      red=0;
      green=0;
      blue=0;
      store.clear();
      closest.clear();
      softShadow.clear();
      
      //Triangle intersection
      for (int i=0;i<num_triangles;i++)
      {
		  triangle = triangles[i];
		  Point A(triangle.v[0].position[0],-triangle.v[0].position[1],triangle.v[0].position[2]);
		  Point B(triangle.v[1].position[0],-triangle.v[1].position[1],triangle.v[1].position[2]);
		  Point C(triangle.v[2].position[0],-triangle.v[2].position[1],triangle.v[2].position[2]);
	  
		  Point n1 = subtract(B,A);
		  Point n2 = subtract(C,A);
		  Point normal = Cross(n1,n2);
		  Normalize(normal);
		  d = Dot(normal,A);
	  
	  	  if (Dot(normal,ray.d) == 0)
	  	  {
	  	  	break;
	  	  }
		  t3 = d/Dot (normal,ray.d);

		  Point Q (ray.d.x*t3,ray.d.y*t3,ray.d.z*t3);
		  
		  //Ordinary method
		  i1 = Dot(Cross(subtract(B,A),subtract(Q,A)),normal);
		  i2 = Dot(Cross(subtract(C,B),subtract(Q,B)),normal);
		  i3 = Dot(Cross(subtract(A,C),subtract(Q,C)),normal);
		  
		  //Barycentric method
		  areaABC = Dot(Cross(n1,n2),normal); 
		  areaPBC = Dot(Cross(subtract(B,Q),subtract(C,Q)),normal);
		  areaPCA = Dot(Cross(subtract(C,Q),subtract(A,Q)),normal);
		 
		  barycentric.x = areaPBC/areaABC;
		  barycentric.y = areaPCA/areaABC;
		  barycentric.z = 1.f - barycentric.x - barycentric.y;
		  if (t3 > 0 &&  barycentric.x > 0 && barycentric.y > 0 && barycentric.z > 0)
		  {
		  	near.intersection = Q;
		  	near.normal = normal;
		  	near.bary = barycentric;
		  	near.triangle = triangle;
		  	closest.push_back(near);
		  }
		}
		if (closest.size() > 0)
		{
			if (closest.size() == 1)
			{
				nearest = closest[0];
			}
			else
			{
				min = closest[0].intersection.z;
				nearest = closest[0];
				for (int i=1;i<closest.size();i++)
				{
					if (closest[i].intersection.z > min)
					{
						min = closest[i].intersection.z;
						nearest = closest[i];
					}
				}
			}
		}
		if (closest.size() > 0)
		{
			Point Q = nearest.intersection;
			barycentric = nearest.bary;
			triangle = nearest.triangle;
		  	shadows = 0;
		  	t4=0;
		  	t5=0;
		    triangleB = true;
		    for (int f=0;f<num_lights;f++)
		    {
		    	shadow = false;
				Light lightSource = lights[f]; 
				Point surface = Q;
				Point light(lightSource.position[0],-lightSource.position[1],lightSource.position[2]);
				shadowRay.p = surface;
				shadowRay.d = subtract(light,surface);
				Normalize (shadowRay.d);
				t5 = (light.x-Q.x)/shadowRay.d.x;
				//Shadow triangle
				for (int x=0;x<num_triangles;x++)
				{
				  tri = triangles[x];
				  Point _A(tri.v[0].position[0],-tri.v[0].position[1],tri.v[0].position[2]);
				  Point _B(tri.v[1].position[0],-tri.v[1].position[1],tri.v[1].position[2]);
				  Point _C(tri.v[2].position[0],-tri.v[2].position[1],tri.v[2].position[2]);
	  
				  Point _n1 = subtract(_B,_A);
				  Point _n2 = subtract(_C,_A);
				  Point _normal = Cross(_n1,_n2);
				  Normalize(_normal);
				  _d = Dot(_normal,_A);
	  
				  if (Dot(_normal,shadowRay.d) == 0)
				  {
					break;
				  }
				  t4 = (_d-Dot(_normal,shadowRay.p))/Dot (_normal,shadowRay.d);
				  Point _Q (shadowRay.p.x+shadowRay.d.x*t4,shadowRay.p.y+shadowRay.d.y*t4,shadowRay.p.z+shadowRay.d.z*t4);
				  _areaABC = Dot(Cross(_n1,_n2),_normal); 
				  _areaPBC = Dot(Cross(subtract(_B,_Q),subtract(_C,_Q)),_normal);
				  _areaPCA = Dot(Cross(subtract(_C,_Q),subtract(_A,_Q)),_normal);
		 
				  bary.x = _areaPBC/_areaABC;
				  bary.y = _areaPCA/_areaABC;
				  bary.z = 1.f - bary.x - bary.y;
				  if (t4 > 0.00001 && t4 < t5 && bary.x > 0 && bary.y > 0 && bary.z > 0)
				  {
					shadow = true;
				  }
				}
				//Shadow sphere
				for (int j=0;j<num_spheres;j++)
				{
					_sphere = spheres[j];
					_b = 2*(shadowRay.d.x*(shadowRay.p.x-_sphere.position[0])+shadowRay.d.y*(shadowRay.p.y-_sphere.position[1]) + shadowRay.d.z*(shadowRay.p.z-_sphere.position[2]));
					_c = pow(shadowRay.p.x-_sphere.position[0],2)+pow(shadowRay.p.y-_sphere.position[1],2)+pow(shadowRay.p.z-_sphere.position[2],2)-_sphere.radius*_sphere.radius;
					_t0 = (-_b+sqrt(_b*_b-4*_c))/2;
					_t1 = (-_b-sqrt(_b*_b-4*_c))/2;
					if ((int)_t0 > 0 || (int)_t1 > 0)
					{
						shadow = true;
					}	
				}
				if (!shadow)
				{
					Point l = subtract(light,surface);	
					if (softshadow)
					{
						for (int s=0;s<100;s++)
						{
							uS = ((double) rand() / (RAND_MAX));
							vS = ((double) rand() / (RAND_MAX));
							qS = 2*(pow(uS-light.x,2)+pow(uS-light.y,2)+pow(uS-light.z,2)-1);
							fS = acos(2*vS-1);
							l = subtract(Point(cos(qS)*sin(fS),sin(qS)*sin(fS),cos(fS)),surface);
							softShadow.push_back(l);
						}
						l.x = 0;
						l.y = 0;
						l.z = 0;
						for (int s=0;s<100;s++)
						{
							l.x += softShadow[s].x;
							l.y += softShadow[s].x;
							l.z += softShadow[s].x;
						}
						l.x /= 100;l.y/=100;l.z/=100;
					}
					Point v(-surface.x,-surface.y,-surface.z);
					nx = barycentric.x*triangle.v[0].normal[0]+barycentric.y*triangle.v[1].normal[0]+barycentric.z*triangle.v[2].normal[0];
					ny = barycentric.x*triangle.v[0].normal[1]+barycentric.y*triangle.v[1].normal[1]+barycentric.z*triangle.v[2].normal[1];
					nz = barycentric.x*triangle.v[0].normal[2]+barycentric.y*triangle.v[1].normal[2]+barycentric.z*triangle.v[2].normal[2];
					Point n(nx,-ny,nz);
					Normalize(l);Normalize(v);Normalize(n);
					scalar = 2*Dot(l,n);
					Point r = subtract(Point(scalar*n.x,scalar*n.y,scalar*n.z),l);
					Normalize(r);
	  
					diffuse = Dot(l,n);
					if (fabs(diffuse) < 1)
					{
						diffuse = fabs(diffuse);
					}
					Shininess = barycentric.x*triangle.v[0].shininess+barycentric.y*triangle.v[1].shininess+barycentric.z*triangle.v[2].shininess;
					spec = Dot(r,v);
					if (fabs(spec) < 1)
					{
						spec = fabs(spec);
					}
					Clamp(diffuse);Clamp(spec);
					specular = pow(spec,Shininess);
			
					redDiffuse = barycentric.x*triangle.v[0].color_diffuse[0]+barycentric.y*triangle.v[1].color_diffuse[0]+barycentric.z*triangle.v[2].color_diffuse[0];
					redSpecular = barycentric.x*triangle.v[0].color_specular[0]+barycentric.y*triangle.v[1].color_specular[0]+barycentric.z*triangle.v[2].color_specular[0];
					greenDiffuse = barycentric.x*triangle.v[0].color_diffuse[1]+barycentric.y*triangle.v[1].color_diffuse[1]+barycentric.z*triangle.v[2].color_diffuse[1];
					greenSpecular = barycentric.x*triangle.v[0].color_specular[1]+barycentric.y*triangle.v[1].color_specular[1]+barycentric.z*triangle.v[2].color_specular[1];
					blueDiffuse = barycentric.x*triangle.v[0].color_diffuse[2]+barycentric.y*triangle.v[1].color_diffuse[2]+barycentric.z*triangle.v[2].color_diffuse[2];
					blueSpecular = barycentric.x*triangle.v[0].color_specular[2]+barycentric.y*triangle.v[1].color_specular[2]+barycentric.z*triangle.v[2].color_specular[2];
			
					red += lightSource.color[0]*(redDiffuse*diffuse+redSpecular*specular);
					green += lightSource.color[1]*(greenDiffuse*diffuse+greenSpecular*specular);
					blue += lightSource.color[2]*(blueDiffuse*diffuse+blueSpecular*specular);
				}
				else
				{
					shadows++;
				}
			  }
			  red += ambient_light[0];
			  green += ambient_light[1];
			  blue += ambient_light[2];
			  Clamp(red);Clamp(green);Clamp(blue);
			  
			  near.z = Q.z;
			  near.color = Point(red,green,blue);
			  store.push_back(near);
	  }
	  
	  red=0;
	  green=0;
	  blue=0;
      //Sphere intersection
      for (int i=0;i<num_spheres;i++)
      {
		  sphere = spheres[i];
		  b = 2*(ray.d.x*(-sphere.position[0])+ray.d.y*(-sphere.position[1]) + ray.d.z*(-sphere.position[2]));
		  c = sphere.position[0]*sphere.position[0]+sphere.position[1]*sphere.position[1]+sphere.position[2]*sphere.position[2]-sphere.radius*sphere.radius;
		  t0 = (-b+sqrt(b*b-4*c))/2;
		  t1 = (-b-sqrt(b*b-4*c))/2;
		  t2 = t1;
		  if (t0 < t1)
		  {
			t2 = t0;
		  }
		  if (t2 > 0)
		  {
			sphereB = true;
			Point surface(ray.d.x * t2,ray.d.y * t2,ray.d.z * t2);
			for (int f=0;f<num_lights;f++)
		    {
		    	shadow = false;
				Light lightSource = lights[f];
				Point light(lightSource.position[0],-lightSource.position[1],lightSource.position[2]);
				shadowRay.p = surface;
				shadowRay.d = subtract(light,surface);
				Normalize (shadowRay.d);
				for (int x=0;x<num_triangles;x++)
				{
				  tri = triangles[x];
				  Point _A(tri.v[0].position[0],-tri.v[0].position[1],tri.v[0].position[2]);
				  Point _B(tri.v[1].position[0],-tri.v[1].position[1],tri.v[1].position[2]);
				  Point _C(tri.v[2].position[0],-tri.v[2].position[1],tri.v[2].position[2]);
	  
				  Point _n1 = subtract(_B,_A);
				  Point _n2 = subtract(_C,_A);
				  Point _normal = Cross(_n1,_n2);
				  Normalize(_normal);
				  _d = Dot(_normal,_A);
	  
				  if (Dot(_normal,shadowRay.d) == 0)
				  {
					break;
				  }
				  t4 = (_d-Dot(_normal,shadowRay.p))/Dot (_normal,shadowRay.d);
				  Point _Q (shadowRay.p.x+shadowRay.d.x*t4,shadowRay.p.y+shadowRay.d.y*t4,shadowRay.p.z+shadowRay.d.z*t4);
				  _areaABC = Dot(Cross(_n1,_n2),_normal); 
				  _areaPBC = Dot(Cross(subtract(_B,_Q),subtract(_C,_Q)),_normal);
				  _areaPCA = Dot(Cross(subtract(_C,_Q),subtract(_A,_Q)),_normal);
		 
				  bary.x = _areaPBC/_areaABC;
				  bary.y = _areaPCA/_areaABC;
				  bary.z = 1.f - bary.x - bary.y;
				  //cout << bary.x << " " << bary.y << " " << bary.z << endl;
				  if (t4 > 0 &&  bary.x > 0 && bary.y > 0 && bary.z > 0)
				  {
					shadow = true;
				  }
				}
				for (int j=0;j<num_spheres;j++)
				{
					_sphere = spheres[j];
					_b = 2*(shadowRay.d.x*(shadowRay.p.x-_sphere.position[0])+shadowRay.d.y*(shadowRay.p.y-_sphere.position[1]) + shadowRay.d.z*(shadowRay.p.z-_sphere.position[2]));
					_c = pow(shadowRay.p.x-_sphere.position[0],2)+pow(shadowRay.p.y-_sphere.position[1],2)+pow(shadowRay.p.z-_sphere.position[2],2)-_sphere.radius*_sphere.radius;
					_t0 = (-_b+sqrt(_b*_b-4*_c))/2;
					_t1 = (-_b-sqrt(_b*_b-4*_c))/2;
					if ((int)_t0 > 0 || (int)_t1 > 0)
					{
						near.z = shadowRay.p.z;
						near.color = Point(0.1,0.1,0.1);
						store.push_back(near);
					}	
				}
				if (!shadow)
				{
					Point l = subtract(light,surface);	
					if (softshadow)
					{
						for (int s=0;s<100;s++)
						{
							uS = ((double) rand() / (RAND_MAX));
							vS = ((double) rand() / (RAND_MAX));
							qS = 2*(pow(uS-light.x,2)+pow(uS-light.y,2)+pow(uS-light.z,2)-1);
							fS = acos(2*vS-1);
							l = subtract(Point(cos(qS)*sin(fS),sin(qS)*sin(fS),cos(fS)),surface);
							softShadow.push_back(l);
						}
						l.x = 0;
						l.y = 0;
						l.z = 0;
						for (int s=0;s<100;s++)
						{
							l.x += softShadow[s].x;
							l.y += softShadow[s].x;
							l.z += softShadow[s].x;
						}
						l.x /= 100;l.y/=100;l.z/=100;
					}
					Point v(-surface.x,-surface.y,-surface.z);
					Point n ((surface.x-sphere.position[0])/sphere.radius,(surface.y-sphere.position[1])/sphere.radius,(surface.z-sphere.position[2])/sphere.radius);
					Normalize(l);Normalize(v);Normalize(n);
					scalar = 2*Dot(l,n);
					Point r = subtract(Point(scalar*n.x,scalar*n.y,scalar*n.z),l);
					Normalize(r);
	  
					diffuse = Dot(l,n);
			
					specular = pow(Dot(r,v),sphere.shininess);
		
					Clamp(diffuse);Clamp(specular);
		
					red += lightSource.color[0]*(sphere.color_diffuse[0]*diffuse+sphere.color_specular[0]*specular);
					green += lightSource.color[1]*(sphere.color_diffuse[1]*diffuse+sphere.color_specular[1]*specular);
					blue += lightSource.color[2]*(sphere.color_diffuse[2]*diffuse+sphere.color_specular[2]*specular);
				}
			}
			red += ambient_light[0];
			green += ambient_light[1];
		    blue += ambient_light[2];
			Clamp(red);Clamp(green);Clamp(blue);
			zSphere = surface.z;
			colorSphere = Point(red,green,blue);
			near.z = zSphere;
			near.color = colorSphere;
			store.push_back(near);
		  }
      }
      if (store.size() > 0)
	  {
	  	if (store.size() == 1)
	  	{
	  		glColor3f(store[0].color.x,store[0].color.y,store[0].color.z);
  	 	    glVertex2i(x,y);
	  	}
	  	else
	  	{
	  		min = store[0].z;
	  		picked = store[0].color;
	  		for (int z=1;z<store.size();z++)
	  		{
	  			if (store[z].z > min)
	  			{
	  				min = store[z].z;
	  				picked = store[z].color;
	  			}
	  		}
	  		glColor3f(picked.x,picked.y,picked.z);
  	 	    glVertex2i(x,y);
	  	}
	  }
      if (!triangleB && !sphereB)
	  {
		  glColor3f(1,1,1);
		  glVertex2i(x,y);
	  }
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

void parse_check(char *expected,char *found)
{
  if(strcasecmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(strcasecmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(strcasecmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(strcasecmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void display()
{

}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();
      if(mode == MODE_JPEG)
	save_jpg();
    }
  once=1;
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
      if (strcmp(argv[2], "softshadow.jpg") == 0)
	  {
		softshadow = true;
	  }
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}