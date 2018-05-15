
#include"fem.h"

# ifndef NOPOISSONCREATE
double jacobien(double x1,double x2,double x3,double y1,double y2, double y3)
{
	return ((x1 - x3)*(y2 - y1) - (x1 - x2)*(y3 - y1));
}

void setVoisin(femPoissonProblem* theProblem)
{
	int a, b, c, temp, flag;
	int i, j, k;
	int * tab = (int *)malloc(sizeof(int) * theProblem->mesh->nElem * 3);
	int loc = 0;

	for (i = 0;i < theProblem->mesh->nElem; i++) {
		a = theProblem->mesh->elem[3 * i];
		b = theProblem->mesh->elem[(3 * i) + 1];
		c = theProblem->mesh->elem[(3 * i) + 2];
		flag = 0;

		for (j = 0; i < (theProblem->mesh->nElem) && flag<3; j++) {
			if (j != i) {
				for (k = 0; k < 3; k++) {
					if (theProblem->mesh->elem[(j * 3) + k] == a) {
						if (theProblem->mesh->elem[j * 3] == b || theProblem->mesh->elem[j * 3] == c || theProblem->mesh->elem[(j * 3) + 1] == b ||
							theProblem->mesh->elem[(j * 3) + 1] == c || theProblem->mesh->elem[(j * 3) + 2] == b || theProblem->mesh->elem[(j * 3) + 2] == c) {
							*(tab + loc) = j;
							loc++;
							flag++;
						}
					}

					else if (theProblem->mesh->elem[(j * 3) + k] == b) {
						if (theProblem->mesh->elem[j * 3] == c || theProblem->mesh->elem[(j * 3) + 1] == c || theProblem->mesh->elem[(j * 3) + 2] == c) {
							*(tab + loc) = j;
							loc++;
							flag++;
						}
					}
				}

				if (j == (theProblem->mesh->nElem - 1) && flag<2) {
					*(tab + loc) = -1;
					flag++;
					loc++;
				}
			}
		}
	}
	theProblem->mesh->voisin = tab;
}

int femInBowl(femPoissonProblem* theProblem, int numElem, double xBowl, double yBowl)
{
	if (numElem == -1)
		return -1;

	double s1 = 0;
	double s2 = 0;
	double s3 = 0;
	double sEq = 0;
	int node[3];

	for (int i = 0; i < 3; i++)
		node[i] = theProblem->mesh->elem[3 * numElem + i];

	double x[] = { theProblem->mesh->X[node[0]], theProblem->mesh->X[node[1]],theProblem->mesh->X[node[2]] };
	double y[] = { theProblem->mesh->Y[node[0]], theProblem->mesh->Y[node[1]],theProblem->mesh->Y[node[2]] };

	s1 = jacobien(xBowl, x[0], x[1], yBowl, y[0], y[1]); 
	s2 = jacobien(xBowl, x[1], x[2], yBowl, y[1], y[2]); 
	s3 = jacobien(xBowl, x[2], x[0], yBowl, y[2], y[0]); 
	return ((s1 >= 0 && s2 >= 0 && s3 >= 0) || (s1 <= 0 && s2 <= 0 && s3 <= 0));
}

//Changer
void isomorphisme(double X[3], double Y[3], double x[2], double *xsi)
{
	double X1 = X[0], X2 = X[1], X3 = X[2];
	double Y1 = Y[0], Y2 = Y[1], Y3 = Y[2];

	xsi[0] = -(X1*Y3 - X3 * Y1 - X1 * x[1] + Y1 * x[0] + X3 * x[1] - Y3 * x[0]) / (X1*Y2 - X2 * Y1 - X1 * Y3 + X3 * Y1 + X2 * Y3 - X3 * Y2);
	xsi[1] = (X1*Y2 - X2 * Y1 - X1 * x[1] + Y1 * x[0] + X2 * x[1] - Y2 * x[0]) / (X1*Y2 - X2 * Y1 - X1 * Y3 + X3 * Y1 + X2 * Y3 - X3 * Y2);
}

int findTriangle(femPoissonProblem* theProblem, double x, double y)
{
	for (int i = 0; i < theProblem->mesh->nElem; i++) {
		if (femInBowl(theProblem, i, x, y) == 1)
			return i;
	}
	return 0;
}

int updateTriangle(femPoissonProblem* theProblem, int iElem, double x, double y, int iter)
{
	int i, j, k;
	int temp;

	if (femInBowl(theProblem, iElem, x, y))
		return iElem;
	else {
		for (int i = 0; i < 3; i++) {
			if (femInBowl(theProblem, theProblem->mesh->voisin[iElem * 3 + i], x, y))
				return (theProblem->mesh->voisin[iElem * 3 + i]);
			temp = theProblem->mesh->voisin[iElem * 3 + i];
			for (int j = 0; i < 3; i++) {
				if (femInBowl(theProblem, theProblem->mesh->voisin[temp * 3 + j], x, y))
					return (theProblem->mesh->voisin[temp * 3 + i]);
			}
		}
		return (findTriangle(theProblem, x, y));
	}
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

void femFullSystemReinit(femPoissonProblem *theProblem)
{
	int size = theProblem->systemX->size;
	for (int i = 0; i < size*(size + 1); i++) {
		theProblem->systemX->B[i] = 0;
		theProblem->systemY->B[i] = 0;
	}
	for (int k = 0; k < size; k++) {
		for (int l = 0; l < size; l++) {
			theProblem->systemX->A[k][l] = 0;
			theProblem->systemY->A[k][l] = 0;
		}
	}
}


# endif
# ifndef NOMESHLOCAL


void femMeshLocal(const femMesh *theMesh, const int i, int *map, double *x, double *y)
{
	int n = theMesh->nLocalNode;
	for (int j = 0; j < n; j++) {
		map[j] = theMesh->elem[(n*i) + j];
		x[j] = theMesh->X[map[j]];
		y[j] = theMesh->Y[map[j]];
	}
}

# endif
# ifndef NOPOISSONSOLVE


void femPoissonSolve(femPoissonProblem *theProblem, femGrains *theGrains)
{
	int localNode = theProblem->mesh->nLocalNode;
	int node = theProblem->mesh->nNode;
	int elem = theProblem->mesh->nElem;
	int map[4];
	double x[4], y[4], phi[4], dphidxsi[4], dphideta[4];

	for (int k = 0; k < elem; k++) {
		femMeshLocal(theProblem->mesh, k, map, x, y);
		for (int p = 0; p < theProblem->rule->n; p++) { 
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
				theProblem->systemX->B[map[f]] += 0;
				theProblem->systemY->B[map[f]] += 0;
			}

			//Partie a decommenter pour resolution de Couette
			/*
			double xsi[3], phi[3], point[2], xs[3], ys[3];
			int elem, map[3];
			for (int i = 0; i < theGrains->n; ++i) {
				elem = theGrains->element[i];
				point[0] = theGrains->x[i];
				point[1] = theGrains->y[i];
				femMeshLocal(theProblem->mesh, elem, map, xs, ys);
				isomorphisme(xs, ys, point, xsi);
				theProblem->space->phi2(xsi[0], xsi[1], phi);
				for (int k = 0; k < 3; ++k) {
					for (int l = 0; l < 3; ++l) {
						theProblem->systemX->A[map[k]][map[l]] += phi[k] * phi[l];
						theProblem->systemY->A[map[k]][map[l]] += phi[k] * phi[l];
					}
					theProblem->systemX->B[map[k]] = phi[k] * theGrains->vx[i];
					theProblem->systemY->B[map[k]] = phi[k] * theGrains->vy[i];
				}
			}
			*/
		}
	}

	double x1, y1;
	for (int z = 0; z < theProblem->edges->nEdge; z++) {
		 x1 = theProblem->mesh->X[theProblem->edges->edges[z].node[0]];
		 y1 = theProblem->mesh->Y[theProblem->edges->edges[z].node[0]];

		if (theProblem->edges->edges[z].elem[1] < 0 && sqrt(pow(x1,2.0)+ pow(y1,2.0)>1.25)) {
			double theta = atan2(y1,x1);
			double speedX = sin(theta);
			double speedY = -cos(theta);
			for (int i = 0; i < 2; i++) {
				femFullSystemConstrain(theProblem->systemX, theProblem->edges->edges[z].node[i], speedX);
				femFullSystemConstrain(theProblem->systemY, theProblem->edges->edges[z].node[i], speedY);
			}

		}

		else if (theProblem->edges->edges[z].elem[1] < 0 ) {
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

	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			double L = 0, gamma = 0, vn = 0;
			double pos[2], nT[2], dv;

			pos[0] = x[j] - x[i];
			pos[1] = y[j] - y[i];

			L = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
			gamma = L - r[i] - r[j];

			nT[0] = pos[0] / L;
			nT[1] = pos[1] / L;
			vn = vx[i] * nT[0] + vy[i] * nT[1] - (vx[j] * nT[0] + vy[j] * nT[1]);
			dv = fmax(0.0, vn + dvContacts[nContacts] - gamma / dt) - dvContacts[nContacts];

			if (iter == 0)
				dv = dvContacts[nContacts];
			else
				dvContacts[nContacts] += dv;

			vx[i] += -dv * nT[0] * m[j] / (m[j] + m[i]);
			vy[i] += -dv * nT[1] * m[j] / (m[j] + m[i]);
			vx[j] += dv * nT[0] * m[i] / (m[j] + m[i]);
			vy[j] += dv * nT[1] * m[i] / (m[j] + m[i]);
			zeta = fmax(zeta, fabs(dv));
			nContacts++;

		}
		double L = 0, gammaOUT = 0, vn = 0, gammaIN = 0;
		double dv = 0, dvOut = 0, dvIn = 0;
		double nT[2];

		L = sqrt(x[i] * x[i] + y[i] * y[i]);

		gammaOUT = rOut - L - r[i];
		gammaIN = L - rIn - r[i];

		nT[0] = x[i] / L;
		nT[1] = y[i] / L;

		vn = vx[i] * nT[0] + vy[i] * nT[1];
		dvOut = fmax(0.0, vn + dvBoundary[i] - gammaOUT / dt);
		dvIn = fmax(0.0, -dvBoundary[i] - vn - gammaIN / dt);
		dv = dvOut - dvIn - dvBoundary[i];

		if (iter == 0)
			dv = dvBoundary[i];
		else
			dvBoundary[i] += dv;

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
	double gamma2 = myGrains->gamma;
	double gx = myGrains->gravity[0];
	double gy = myGrains->gravity[1];


	
	int elem, map[3];
	double xloc[3], yloc[3], phi[3], xsi[2], p[2], ux, uy;
	for (int i = 0; i < n; ++i) {
		elem = myGrains->element[i];
		if (elem < 0) {
			ux = 0;
			uy = 0;
		}
		else {
			p[0] = x[i];
			p[1] = y[i];
			femMeshLocal(theProblem->mesh, elem, map, xloc, yloc);
			isomorphisme(xloc, yloc, p, xsi);
			ux = 0;
			uy = 0;
			theProblem->space->phi2(xsi[0], xsi[1], phi);
			for (int j = 0; j < 3; ++j) {
				ux += phi[j] * theProblem->systemX->B[map[j]];
				uy += phi[j] * theProblem->systemY->B[map[j]];
			}
		}
		vx[i] += (m[i] * gx -   (vx[i] - ux)) * dt / m[i];
		vy[i] += (m[i] * gy - (vy[i] - uy)) * dt / m[i];
	}

	do {
		zeta = femGrainsContactIterate(myGrains, dt, iter);
		iter++;
	} while ((zeta > tol / dt && iter < iterMax) || iter == 1);
	printf("iterations = %4d : error = %14.7e \n", iter - 1, zeta);


	for (i = 0; i < n; ++i) {
		x[i] += vx[i] * dt;
		y[i] += vy[i] * dt;
	}
	for (int k = 0; k < n; k++)
		myGrains->element[k] = updateTriangle(theProblem ,myGrains->element[k], x[k], y[k], 0);
}


#endif


# endif