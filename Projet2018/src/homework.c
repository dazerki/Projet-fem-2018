
#include"fem.h"

# ifndef NOPOISSONCREATE

int femInBowl(femPoissonProblem* theProblem, int numElem, int xBowl, int yBowl)
{
	double s = 0;
	double sEq = 0;
	int node[] = { theProblem->mesh->elem[numElem], theProblem->mesh->elem[numElem + 1],theProblem->mesh->elem[numElem + 2] };
	double x[] = { theProblem->mesh->X[node[0]], theProblem->mesh->X[node[1]],theProblem->mesh->X[node[2]] };
	double y[] = { theProblem->mesh->Y[node[0]], theProblem->mesh->Y[node[1]],theProblem->mesh->Y[node[2]] };
	s = 0.5 * fabs((x[1] - xBowl)*(y[2] - yBowl) - (x[2] - xBowl)*(y[1] - yBowl));
	s += 0.5 * fabs((x[0] - xBowl)*(y[1] - yBowl) - (x[1] - xBowl)*(y[0] - yBowl));
	s += 0.5 * fabs((x[2] - xBowl)*(y[0] - yBowl) - (x[0] - xBowl)*(y[2] - yBowl));
	sEq = 0.5 * fabs((x[1] - x[0])*(y[2] - y[0]) - (x[2] - x[0])*(y[1] - y[0]));
	if (s<=sEq)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
int findTriangle(femPoissonProblem* theProblem, double x, double y) {
	for (int i = 0; i < theProblem->mesh->nElem; i++) {
		if (femInBowl(theProblem, i, x, y) == 1)
			return i;
	}
	return 0;
}

femPoissonProblem *femPoissonCreate(const char *filename)
{
	femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
	theProblem->mesh = femMeshRead(filename);
	theProblem->edges = femEdgesCreate(theProblem->mesh);
	theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
	theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
	
	theProblem->systemX = femFullSystemCreate(theProblem->mesh->nNode);
	theProblem->systemY = femFullSystemCreate(theProblem->mesh->nNode);
	return theProblem;
}

# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
	femFullSystemFree(theProblem->systemX);
	femFullSystemFree(theProblem->systemY);
	femIntegrationFree(theProblem->rule);
	femDiscreteFree(theProblem->space);
	femEdgesFree(theProblem->edges);
	femMeshFree(theProblem->mesh);
	free(theProblem);
}


# endif
# ifndef NOMESHLOCAL


void femMeshLocal(const femMesh *theMesh, const int i, int *map, double *x, double *y)
{
	int n = theMesh->nLocalNode;
	for (int j = 0; j < n; j++)
	{
		map[j] = theMesh->elem[(n*i) + j];
		x[j] = theMesh->X[map[j]];
		y[j] = theMesh->Y[map[j]];
	}
}

# endif
# ifndef NOPOISSONSOLVE


void femPoissonSolve(femPoissonProblem *theProblem)
{

	// Calcul de A: la matrice de raideur
	// Aij=1/2 \Sum{j}( \Sum{i}( Ui Uj <\Nabla \tau_i, \Nabla \tau_j >))
	// C'est parti pour l'amusement -_-
	int localNode = theProblem->mesh->nLocalNode;
	int node = theProblem->mesh->nNode;
	int elem = theProblem->mesh->nElem;
	int map[4];
	double x[4], y[4], phi[4], dphidxsi[4], dphideta[4];
	for (int k = 0; k < elem; k++) {
		femMeshLocal(theProblem->mesh, k, map, x, y);
		for (int p = 0; p < theProblem->rule->n; p++) { // rule->n vaut 4? 
														//theProblem->space->x2(xsi, eta);
			theProblem->space->phi2(theProblem->rule->xsi[p], theProblem->rule->eta[p], phi);
			theProblem->space->dphi2dx(theProblem->rule->xsi[p], theProblem->rule->eta[p], dphidxsi, dphideta);
			double dxdxsi = 0, dxdeta = 0, dydxsi = 0, dydeta = 0;
			for (int w = 0; w < theProblem->space->n; w++) {

				dxdxsi = dxdxsi + x[w] * dphidxsi[w];
				dydxsi = dydxsi + y[w] * dphidxsi[w];
				dxdeta = dxdeta + x[w] * dphideta[w];
				dydeta = dydeta + y[w] * dphideta[w];
			}
			double Je = fabs(dxdxsi * dydeta - dydxsi * dxdeta);

			double dphidx[4], dphidy[4];
			for (int i = 0; i < theProblem->space->n; i++) {
				dphidx[i] = (dydeta * dphidxsi[i] - dydxsi * dphideta[i]) / Je;
				dphidy[i] = (dxdxsi * dphideta[i] - dxdeta * dphidxsi[i]) / Je;
			}
			for (int i = 0; i < theProblem->space->n; i++) {
				for (int j = 0; j < theProblem->space->n; j++) {
					theProblem->systemX->A[map[i]][map[j]] +=  (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * Je * theProblem->rule->weight[p];
					theProblem->systemY->A[map[i]][map[j]] +=  (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * Je * theProblem->rule->weight[p];
				}
			}
			for (int f = 0; f < theProblem->space->n; f++) {
				theProblem->systemX->B[map[f]] = 0; //theProblem->system->Bx[map[f]] + phi[f] * Je *  theProblem->rule->weight[p];
				theProblem->systemY->B[map[f]] = 0;
			}
			
		}
	}
	double x1, y1;
	for (int z = 0; z < theProblem->edges->nEdge; z++) {
		 x1 = theProblem->mesh->X[theProblem->edges->edges[z].node[0]];
		 y1 = theProblem->mesh->Y[theProblem->edges->edges[z].node[0]];
		if (theProblem->edges->edges[z].elem[1] < 0 && sqrt(pow(x1,2.0)+ pow(y1,2.0)>1.25)){
			double theta = atan2(y1,x1);
			double speedX = sin(theta);
			double speedY = -cos(theta);
			printf("%f, %f, %f, %f \n",x1, y1, speedX, speedY);
			for (int i = 0; i < 2; i++) {
				femFullSystemConstrain(theProblem->systemX, theProblem->edges->edges[z].node[i], speedX);
				femFullSystemConstrain(theProblem->systemY, theProblem->edges->edges[z].node[i], speedY);
			}

		}
		else if (theProblem->edges->edges[z].elem[1] < 0 )
		{
			for (int i = 0; i < 2; i++) {
				femFullSystemConstrain(theProblem->systemX, theProblem->edges->edges[z].node[i], 0.0);
				femFullSystemConstrain(theProblem->systemY, theProblem->edges->edges[z].node[i], 0.0);
			}
		}
		

	}
	femFullSystemEliminate(theProblem->systemX);
	femFullSystemEliminate(theProblem->systemY);
	
}


#ifndef NOCONTACTITERATE

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)
{
	int n = myGrains->n;
	double *x = myGrains->x;
	double *y = myGrains->y;
	double *m = myGrains->m;
	double *r = myGrains->r;
	double *vy = myGrains->vy;
	double *vx = myGrains->vx;
	double *dvBoundary = myGrains->dvBoundary;
	double *dvContacts = myGrains->dvContacts;
	double rIn = myGrains->radiusIn;
	double rOut = myGrains->radiusOut;

	double zeta = 0.0;
	int nContacts = 0;

	//
	//  A FAIRE.... :-)    Difficile, difficile :-)  OOOH OUI!
	//
	/*int k = 0;
	for (int i = 0; i < n; i++) {
	dvBoundary[i] = 0;
	for (int j = i+1; j < n; j++)
	{
	dvContacts[k] = 0;
	k++;
	}
	}*/
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			double L = 0, gamma = 0, vn = 0;
			double pos[2], nT[2];
			pos[0] = x[j] - x[i];
			pos[1] = y[j] - y[i];
			L = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
			gamma = L - r[i] - r[j];
			nT[0] = pos[0] / L;
			nT[1] = pos[1] / L;
			vn = vx[i] * nT[0] + vy[i] * nT[1] - (vx[j] * nT[0] + vy[j] * nT[1]);
			double dv;
			dv = fmax(0.0, vn + dvContacts[nContacts] - gamma / dt) - dvContacts[nContacts];
			if (iter == 0) {
				dv = dvContacts[nContacts];
			}
			else
			{
				dvContacts[nContacts] += dv;
			}
			vx[i] += -dv * nT[0] * m[j] / (m[j] + m[i]);
			vy[i] += -dv * nT[1] * m[j] / (m[j] + m[i]);
			vx[j] += dv * nT[0] * m[i] / (m[j] + m[i]);
			vy[j] += dv * nT[1] * m[i] / (m[j] + m[i]);
			zeta = fmax(zeta, fabs(dv));
			nContacts++;

		}
		double L = 0, gammaOUT = 0, vn = 0, gammaIN = 0;
		double nT[2];
		L = sqrt(x[i] * x[i] + y[i] * y[i]);
		gammaOUT = rOut - L - r[i];
		gammaIN = L - rIn - r[i];
		nT[0] = x[i] / L;
		nT[1] = y[i] / L;
		vn = vx[i] * nT[0] + vy[i] * nT[1];
		double dv = 0, dvOut = 0, dvIn = 0;
		dvOut = fmax(0.0, vn + dvBoundary[i] - gammaOUT / dt);
		dvIn = fmax(0.0, -dvBoundary[i] - vn - gammaIN / dt);
		dv = dvOut - dvIn - dvBoundary[i];
		/*if (-dvBoundary[i] - vn - gammaIN / dt > 0) {
		dv = -di - dvBoundary[i];
		}
		else if (vn + dvBoundary[i] - gammaOUT / dt > 0) {
		dv = di - dvBoundary[i];
		}
		else {
		dv = 0;
		}*/
		if (iter == 0) {
			dv = dvBoundary[i];
		}
		else
		{
			dvBoundary[i] += dv;
		}
		vx[i] += -dv * nT[0];
		vy[i] += -dv * nT[1];
		zeta = fmax(zeta, fabs(dv));
	}


	return zeta;

}

#endif
#ifndef NOUPDATE

void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax, femPoissonProblem *theProblem)
{
	int n = myGrains->n;
	int i, iter = 0;
	double zeta;
	double *x = myGrains->x;
	double *y = myGrains->y;
	double *m = myGrains->m;
	double *vy = myGrains->vy;
	double *vx = myGrains->vx;
	double gamma = myGrains->gamma;
	double gx = myGrains->gravity[0];
	double gy = myGrains->gravity[1];


	// 
	// -1- Calcul des nouvelles vitesses des grains sur base de la gravité et de la trainee
	//
	for (int j = 0; j < n; j++) {
		int elem = findTriangle(theProblem, x[j], y[j]);
		vx[j] +=  - (dt * gamma/m[j])*(vx[j]- theProblem->systemX->B[elem] );
		vy[j] += dt * gy - (dt * gamma/m[j])*(vy[j]- theProblem->systemY->B[elem]) ;
	}
	printf("%f", gamma);
	//
	// -2- Correction des vitesses pour tenir compte des contacts        
	//       
	do {
		zeta = femGrainsContactIterate(myGrains, dt, iter);
		iter++;
	} while ((zeta > tol / dt && iter < iterMax) || iter == 1);
	printf("iterations = %4d : error = %14.7e \n", iter - 1, zeta);

	//  
	// -3- Calcul des nouvelles positions sans penetrations de points entre eux
	//
	for (i = 0; i < n; ++i) {
		x[i] += vx[i] * dt;
		y[i] += vy[i] * dt;
	}
}


#endif


# endif
