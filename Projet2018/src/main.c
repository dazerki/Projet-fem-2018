/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

#include "glfem.h"



int main(void)
{   
 
    femPoissonProblem* theProblem = femPoissonCreate("C:\\Users\\antoi\\Desktop\\UCL\\Q4\\myFem-Poisson\\data\\projet_medium.txt");
    
    // Pour Windows, remplacer l'argument :
    // ("../data/triangles_166.txt") 
    // par :
    // ("..\\data\\triangles_166.txt") 
    //
    // Sorry for the inconvenience :-)
    // On r�fl�chit pour rendre cela plus transparent dans les homeworks suivants :-)
    // Be patient !
    
    
    printf("Number of elements    : %4d\n", theProblem->mesh->nElem);
    printf("Number of local nodes : %4d\n", theProblem->mesh->nLocalNode);
    printf("Number of segments    : %4d\n", theProblem->edges->nBoundary);
    printf("Number of unknowns    : %4d\n", theProblem->systemX->size);

	//int l = femInBowl(theProblem, 0, 4, 5);
	//printf("%d", l);

    femPoissonSolve(theProblem);   
	int    n = 15;
	double radius = 0.1;
	double mass = 0.1;
	double radiusIn = 0.5;
	double radiusOut = 2.0;
	double dt = 1.0e-1;
	double tEnd = 8.0;
	double tol = 1e-6;
	double t = 0;
	double iterMax = 100;
	femGrains* theGrains = femGrainsCreateSimple(n, radius, mass, radiusIn, radiusOut);
	int quelTriangle = femInBowl(theProblem,0, 1.50283, 0.7032);
	printf("%d \n", quelTriangle);
	//glfwMakeContextCurrent(window);
	double *V = malloc(sizeof(double) * theProblem->systemX->size);
	for (int i = 0; i < theProblem->systemX->size; i++) {
		V[i] = sqrt(pow(theProblem->systemX->B[i], 2) + pow(theProblem->systemY->B[i], 2));
	}
	int theRunningMode = 1.0;
	float theVelocityFactor = 0.25;
 
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
		glfemReshapeWindows2(radiusOut, w, h);
        glfemPlotField(theProblem->mesh,V);            
		
		
		for (i = 0; i < theGrains->n; i++) {

			glColor3f(1, 0, 0);
			glfemDrawDisk(theGrains->x[i], theGrains->y[i], theGrains->r[i]);
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
				_sleep(100);
				femGrainsUpdate(theGrains, dt, tol, iterMax,theProblem);
				t += dt;

			}
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
