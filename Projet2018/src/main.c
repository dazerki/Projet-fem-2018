/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

#include <time.h>
#include "glfem.h"



int main(void)
{   
 
	femPoissonProblem* theProblem = femPoissonCreate("C:\\Users\\Louis\\Desktop\\Projet2018\\data\\projet_medium.txt");
    //femPoissonProblem* theProblem = femPoissonCreate("C:\\Users\\antoi\\Desktop\\UCL\\Q4\\myFem-Poisson\\data\\projet_medium.txt");
    
    // Pour Windows, remplacer l'argument :
    // ("../data/triangles_166.txt") 
    // par :
    // ("..\\data\\triangles_166.txt") 
    //
    // Sorry for the inconvenience :-)
    // On réfléchit pour rendre cela plus transparent dans les homeworks suivants :-)
    // Be patient !
    
    
    printf("Number of elements    : %4d\n", theProblem->mesh->nElem);
    printf("Number of local nodes : %4d\n", theProblem->mesh->nLocalNode);
    printf("Number of segments    : %4d\n", theProblem->edges->nBoundary);
    printf("Number of unknowns    : %4d\n", theProblem->systemX->size);

	//int l = femInBowl(theProblem, 0, 4, 5);
	//printf("%d", l);

      

	int    n = 15;
	double radius = 0.05;
	double mass = 0.05;
	double radiusIn = 0.4;
	double radiusOut = 2.0;
	double dt = 1.0e-1;
	double tEnd = 20.0;
	double tol = 1e-6;
	double t = 0;
	double iterMax = 100;
	femGrains* theGrains = femGrainsCreateSimple(n, radius, mass, radiusIn, radiusOut);
	for (int k = 0; k < n; k++) {
		theGrains->element[k] = findTriangle(theProblem, theGrains->x[k], theGrains->y[k]);
	}
	setVoisin(theProblem);
	
	//int quelTriangle = findTriangle(theProblem, -1.65, 0.9);
	//printf("%d \n", quelTriangle);
	//glfwMakeContextCurrent(window);
	double *V = malloc(sizeof(double) * theProblem->systemX->size);
	femPoissonSolve(theProblem, theGrains);
	int theRunningMode = 1.0;
	float theVelocityFactor = 1.0;
 
    printf("Maximum value : %.4f\n", femMax(theProblem->systemX->B,theProblem->systemX->size));
    fflush(stdout);
    
    char theMessage[256];
    sprintf(theMessage, "Max : %.4f", femMax(theProblem->systemX->B,theProblem->systemX->size));
    
    GLFWwindow* window = glfemInit("MECA1120 : Projet ");
    glfwMakeContextCurrent(window);
    do {
		int i, w, h;
		double currentTime = glfwGetTime();

        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theProblem->mesh,w,h);
		//glfemReshapeWindows2(radiusOut, w, h);
		       
		for (int i = 0; i < theProblem->systemX->size; i++) {
			V[i] = sqrt(pow(theProblem->systemX->B[i], 2) + pow(theProblem->systemY->B[i], 2));
		}
		glfemPlotField(theProblem->mesh, V);
		
		for (i = 0; i < theGrains->n; i++) {

			glColor3f(1, 1, 1);
			glfemDrawDisk(theGrains->x[i], theGrains->y[i], theGrains->r[i]);
		}

			glColor3f(0, 0, 0); glfemDrawCircle(0, 0, radiusOut);
			glColor3f(0, 0, 0); glfemDrawCircle(0, 0, radiusIn);
			
			char theMessage[256];
			sprintf(theMessage, "Time = %g sec", t);
			glColor3f(1.0, 0.0, 0.0); glfemDrawMessage(20, 460, theMessage);
			glfwSwapBuffers(window);
			glfwPollEvents();

		
			if (t < tEnd && theRunningMode == 1) {
				printf("Time = %4g : ", t);
				
				// A decommenter pour pouvoir progresser pas par pas
				//       printf("press CR to compute the next time step >>");
				//      char c= getchar();
				//
				//_sleep(100);
				//femFullSystemReinit(theProblem);
				//femPoissonSolve(theProblem, theGrains);
				femGrainsUpdate(theGrains, dt, tol, iterMax,theProblem);
				t += dt;

			}
			
		
			while (glfwGetTime() - currentTime < theVelocityFactor) {
				if (glfwGetKey(window, 'R') == GLFW_PRESS)
					theRunningMode = 1;
				if (glfwGetKey(window, 'S') == GLFW_PRESS)
					theRunningMode = 0;
			}
			
		
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             (!glfwWindowShouldClose(window)));

	
    
	glfwTerminate();
	femGrainsFree(theGrains);
	femPoissonFree(theProblem);
	exit(EXIT_SUCCESS);
    
    

}

