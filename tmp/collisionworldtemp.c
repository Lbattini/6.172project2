/** 
 * collision_world.c -- detect and handle line segment intersections
 * Copyright (c) 2012 the Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE. 
 **/

#include "./collision_world.h"
#include<cilk/cilk.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "./intersection_detection.h"
#include "./intersection_event_list.h"
#include "./line.h"
#define N 10000
#define R 20
int quadtree[N][20*R];
bool isinsqr(Vec p,double xu,double yu, double xl, double yl){
  if(p.x>=xu && p.x<=xl && p.y>=yu && p.y<=yl)
    return 1;
  return 0;
}
void initquadtree(int i,double xu,double yu, double xl,double yl,CollisionWorld* collisionWorld){
  int j,k,child,lim;
  bool moved[10*R];
  double nxu,nyu,nxl,nyl;
  double t = collisionWorld->timeStep;
  //printf("%d %d\n",i,quadtree[i][0]);
  //printf("%f %f %f %f\n",xu,yu,xl,yl);
  if(quadtree[i][0]>=R){
    //startsz=quadtree[i][0];
    //for(j=1;j<=quadtree[i][0];j++) 
      //moved[j]=0;
    lim=quadtree[i][0];
    for(j=1;j<=lim;j++){
      Line *line=collisionWorld->lines[quadtree[i][j]];
      Vec shift,p3,p4;
      moved[j]=0;
      shift=Vec_multiply(line->velocity, t);
      p3= Vec_add(line->p1, shift);
      p4= Vec_add(line->p2, shift);
      nxu=xu;
      nyu=yu;
      nxl=(xu+xl)/2;
      nyl=(yu+yl)/2;
      if(isinsqr(line->p1,nxu,nyu,nxl,nyl) && isinsqr(line->p2,nxu,nyu,nxl,nyl) && isinsqr(p3,nxu,nyu,nxl,nyl) && isinsqr(p4,nxu,nyu,nxl,nyl)){
        child=4*i-2;
        if(quadtree[child][0]==-1)
          quadtree[child][0]=0;
        quadtree[child][0]++;
        quadtree[child][quadtree[child][0]]=quadtree[i][j];
        moved[j]=1;
        if(quadtree[child][0]>=R)
          initquadtree(child,nxu,nyu,nxl,nyl,collisionWorld);//to do for other 3
      }
      else{
      nxu=(xu+xl)/2;
      nyu=yu;
      nxl=xl;
      nyl=(yu+yl)/2;
      if(isinsqr(line->p1,nxu,nyu,nxl,nyl) && isinsqr(line->p2,nxu,nyu,nxl,nyl) && isinsqr(p3,nxu,nyu,nxl,nyl) && isinsqr(p4,nxu,nyu,nxl,nyl)){
        child=4*i-1;
        if(quadtree[child][0]==-1)
          quadtree[child][0]=0;
         quadtree[child][0]++;
        quadtree[child][quadtree[child][0]]=quadtree[i][j];
        moved[j]=1;
        if(quadtree[child][0]>=R)
        initquadtree(child,nxu,nyu,nxl,nyl,collisionWorld);
      }
      else{
      nxu=xu;
      nyu=(yu+yl)/2;
      nxl=(xu+xl)/2;
      nyl=yl;
      if(isinsqr(line->p1,nxu,nyu,nxl,nyl) && isinsqr(line->p2,nxu,nyu,nxl,nyl) && isinsqr(p3,nxu,nyu,nxl,nyl) && isinsqr(p4,nxu,nyu,nxl,nyl)){
        child=4*i;
        if(quadtree[child][0]==-1)
          quadtree[child][0]=0;
        quadtree[child][0]++;
        quadtree[child][quadtree[child][0]]=quadtree[i][j];
        moved[j]=1;
        if(quadtree[child][0]>=R)
        initquadtree(child,nxu,nyu,nxl,nyl,collisionWorld);
      }
      else{
      nxu=(xu+xl)/2;
      nyu=(yu+yl)/2;
      nxl=xl;
      nyl=yl;
      if(isinsqr(line->p1,nxu,nyu,nxl,nyl) && isinsqr(line->p2,nxu,nyu,nxl,nyl) && isinsqr(p3,nxu,nyu,nxl,nyl) && isinsqr(p4,nxu,nyu,nxl,nyl)){
        child=4*i+1;
        //printf("d");
        if(quadtree[child][0]==-1)
          quadtree[child][0]=0;
        quadtree[child][0]++;
        quadtree[child][quadtree[child][0]]=quadtree[i][j];
        moved[j]=1;
        if(quadtree[child][0]>=R){
        initquadtree(child,nxu,nyu,nxl,nyl,collisionWorld);}
      }}}}
    }
    for(j=1,k=1;j<=quadtree[i][0];j++){
      if(!moved[j]){
        quadtree[i][k]=quadtree[i][j];
        k++;
      }
    }
    //printf("%d %d quadtree k \n",quadtree[i][0],k);
    quadtree[i][0]=k-1;
  }
}
void quadintersection(int i,CollisionWorld* collisionWorld,IntersectionEventList* intersectionEventList){
  if(quadtree[i][0]==-1)
    return;
  int j,k,parent;
  IntersectionEventList* intersectioncopy;
  //first,check intersection in the same quadrant, at the same level
  for (j = 1; j <= quadtree[i][0]; j++) {
    Line *l1 = collisionWorld->lines[quadtree[i][j]];
   
    for (k = j + 1; k <= quadtree[i][0]; k++) {
      Line *l2 = collisionWorld->lines[quadtree[i][k]];
      // intersect expects compareLines(l1, l2) < 0 to be true.
      // Swap l1 and l2, if necessary.
      if (compareLines(l1, l2) >= 0) {
        Line *temp = l1;
        l1 = l2;
        l2 = temp;
      }

      IntersectionType intersectionType =
          intersect(l1, l2, collisionWorld->timeStep);
      if (intersectionType != NO_INTERSECTION) {
        IntersectionEventList_appendNode(intersectionEventList, l1, l2,
                                         intersectionType);
        collisionWorld->numLineLineCollisions++;
      }
    }
  }
  //check intersection with bigger quadrants( parents), up to the root
  for (j = 1; j <= quadtree[i][0]; j++) {
    //Line *l1 = collisionWorld->lines[quadtree[i][j]];
    parent=(i+2)/4;
    while(parent>=1){
      for (k = 1; k <= quadtree[parent][0]; k++) {
        Line *l2 = collisionWorld->lines[quadtree[parent][k]];
        Line *l1 = collisionWorld->lines[quadtree[i][j]];
        // intersect expects compareLines(l1, l2) < 0 to be true.
        // Swap l1 and l2, if necessary.
        if (compareLines(l1, l2) >= 0) {
          Line *temp = l1;
          l1 = l2;
          l2 = temp;
        }

        IntersectionType intersectionType =
          intersect(l1, l2, collisionWorld->timeStep);
        if (intersectionType != NO_INTERSECTION) {
          IntersectionEventList_appendNode(intersectionEventList, l1, l2,
                                         intersectionType);
          collisionWorld->numLineLineCollisions++;
        }
      }
      parent=(parent+2)/4;
    }
  }
  intersectioncopy=intersectionEventList;
  quadintersection(4*i-2,collisionWorld,intersectioncopy);
  intersectionEventList=intersectioncopy;
  quadintersection(4*i-1,collisionWorld,intersectionEventList);
  //cilk_sync;
  quadintersection(4*i,collisionWorld,intersectionEventList);
  quadintersection(4*i+1,collisionWorld,intersectionEventList);
}
CollisionWorld* CollisionWorld_new(const unsigned int capacity) {
  assert(capacity > 0);

  CollisionWorld* collisionWorld = malloc(sizeof(CollisionWorld));
  if (collisionWorld == NULL) {
    return NULL;
  }

  collisionWorld->numLineWallCollisions = 0;
  collisionWorld->numLineLineCollisions = 0;
  collisionWorld->timeStep = 0.5;
  collisionWorld->lines = malloc(capacity * sizeof(Line*));
  collisionWorld->numOfLines = 0;
  return collisionWorld;
}

void CollisionWorld_delete(CollisionWorld* collisionWorld) {
  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    free(collisionWorld->lines[i]);
  }
  free(collisionWorld->lines);
  free(collisionWorld);
}

unsigned int CollisionWorld_getNumOfLines(CollisionWorld* collisionWorld) {
  return collisionWorld->numOfLines;
}

void CollisionWorld_addLine(CollisionWorld* collisionWorld, Line *line) {
  collisionWorld->lines[collisionWorld->numOfLines] = line;
  collisionWorld->numOfLines++;
}

Line* CollisionWorld_getLine(CollisionWorld* collisionWorld,
                             const unsigned int index) {
  if (index >= collisionWorld->numOfLines) {
    return NULL;
  }
  return collisionWorld->lines[index];
}

void CollisionWorld_updateLines(CollisionWorld* collisionWorld) {
  CollisionWorld_detectIntersection(collisionWorld);
  CollisionWorld_updatePosition(collisionWorld);
  CollisionWorld_lineWallCollision(collisionWorld);
}

void CollisionWorld_updatePosition(CollisionWorld* collisionWorld) {
  double t = collisionWorld->timeStep;
  Vec shift;
  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    Line *line = collisionWorld->lines[i];
    shift=Vec_multiply(line->velocity, t);
    line->p1 = Vec_add(line->p1, shift);
    line->p2 = Vec_add(line->p2, shift);
  }
}

void CollisionWorld_lineWallCollision(CollisionWorld* collisionWorld) {
  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    Line *line = collisionWorld->lines[i];
    bool collide = false;

    // Right side
    if ((line->p1.x > BOX_XMAX || line->p2.x > BOX_XMAX)
        && (line->velocity.x > 0)) {
      line->velocity.x = -line->velocity.x;
      collide = true;
    }
    // Left side
    if ((line->p1.x < BOX_XMIN || line->p2.x < BOX_XMIN)
        && (line->velocity.x < 0)) {
      line->velocity.x = -line->velocity.x;
      collide = true;
    }
    // Top side
    if ((line->p1.y > BOX_YMAX || line->p2.y > BOX_YMAX)
        && (line->velocity.y > 0)) {
      line->velocity.y = -line->velocity.y;
      collide = true;
    }
    // Bottom side
    if ((line->p1.y < BOX_YMIN || line->p2.y < BOX_YMIN)
        && (line->velocity.y < 0)) {
      line->velocity.y = -line->velocity.y;
      collide = true;
    }
    // Update total number of collisions.
    if (collide == true) {
      collisionWorld->numLineWallCollisions++;
    }
  }
}

void CollisionWorld_detectIntersection(CollisionWorld* collisionWorld) {
  //printf("_____\n");
  IntersectionEventList intersectionEventList = IntersectionEventList_make();
  for (int i = 0; i < N; i++) quadtree[i][0]=-1;
  //initialize quadtree, by adding 3*R lines to the root, and then pushing them by calling initquadtree
  int j=1;
  for (int i = 0; i < collisionWorld->numOfLines; i++,j++){
    quadtree[1][j]=i; 
    if(j-quadtree[1][0]>=R){
      quadtree[1][0]=j;
      initquadtree(1,0.5,0.5,1,1,collisionWorld);
      j=quadtree[1][0];
    }
  }
  quadtree[1][0]=j-1;
  initquadtree(1,0.5,0.5,1,1,collisionWorld);
  quadintersection(1,collisionWorld,&intersectionEventList);
  // Test all line-line pairs to see if they will intersect before the
  // next time step.
  /*
  for (int i = 0; i < collisionWorld->numOfLines; i++) {
    Line *l1 = collisionWorld->lines[i];
    for (int j = i + 1; j < collisionWorld->numOfLines; j++) {
      Line *l2 = collisionWorld->lines[j];
      // intersect expects compareLines(l1, l2) < 0 to be true.
      // Swap l1 and l2, if necessary.
      if (compareLines(l1, l2) >= 0) {
        Line *temp = l1;
        l1 = l2;
        l2 = temp;
      }
      IntersectionType intersectionType =
          intersect(l1, l2, collisionWorld->timeStep);
      if (intersectionType != NO_INTERSECTION) {
        IntersectionEventList_appendNode(&intersectionEventList, l1, l2,
                                         intersectionType);
        collisionWorld->numLineLineCollisions++;
      }
    }
  }
  */
  // Sort the intersection event list.
  IntersectionEventNode* startNode = intersectionEventList.head;
  while (startNode != NULL) {
    IntersectionEventNode* minNode = startNode;
    IntersectionEventNode* curNode = startNode->next;
    while (curNode != NULL) {
      if (IntersectionEventNode_compareData(curNode, minNode) < 0) {
        minNode = curNode;
      }
      curNode = curNode->next;
    }
    if (minNode != startNode) {
      IntersectionEventNode_swapData(minNode, startNode);
    }
    startNode = startNode->next;
  }

  // Call the collision solver for each intersection event.
  IntersectionEventNode* curNode = intersectionEventList.head;

  while (curNode != NULL) {
    //printf("%f %f %f %f\n",curNode->l1->p1.x,curNode->l1->p1.y,curNode->l2->p1.x,curNode->l2->p1.y);
    CollisionWorld_collisionSolver(collisionWorld, curNode->l1, curNode->l2,
                                   curNode->intersectionType);
    curNode = curNode->next;
  }

  IntersectionEventList_deleteNodes(&intersectionEventList);
}

unsigned int CollisionWorld_getNumLineWallCollisions(
    CollisionWorld* collisionWorld) {
  return collisionWorld->numLineWallCollisions;
}

unsigned int CollisionWorld_getNumLineLineCollisions(
    CollisionWorld* collisionWorld) {
  return collisionWorld->numLineLineCollisions;
}

void CollisionWorld_collisionSolver(CollisionWorld* collisionWorld,
                                    Line *l1, Line *l2,
                                    IntersectionType intersectionType) {
  assert(compareLines(l1, l2) < 0);
  assert(intersectionType == L1_WITH_L2
         || intersectionType == L2_WITH_L1
         || intersectionType == ALREADY_INTERSECTED);

  // Despite our efforts to determine whether lines will intersect ahead
  // of time (and to modify their velocities appropriately), our
  // simplified model can sometimes cause lines to intersect.  In such a
  // case, we compute velocities so that the two lines can get unstuck in
  // the fastest possible way, while still conserving momentum and kinetic
  // energy.
  if (intersectionType == ALREADY_INTERSECTED) {
    Vec p = getIntersectionPoint(l1->p1, l1->p2, l2->p1, l2->p2);
    Vec l1p1,l1p2,l2p1,l2p2;
    l1p1=Vec_subtract(l1->p1, p);
    l1p2=Vec_subtract(l1->p2, p);
    l2p1=Vec_subtract(l2->p1, p);
    l2p2=Vec_subtract(l2->p2, p);
    if (Vec_length(l1p1)
        < Vec_length(l1p2)) {
      l1->velocity = Vec_multiply(Vec_normalize(l1p2),
                                  Vec_length(l1->velocity));
    } else {
      l1->velocity = Vec_multiply(Vec_normalize(l1p1),
                                  Vec_length(l1->velocity));
    }
    if (Vec_length(l2p1)
        < Vec_length(l2p2)) {
      l2->velocity = Vec_multiply(Vec_normalize(l2p2),
                                  Vec_length(l2->velocity));
    } else {
      l2->velocity = Vec_multiply(Vec_normalize(l2p1),
                                  Vec_length(l2->velocity));
    }
    return;
  }

  // Compute the collision face/normal vectors.
  Vec face;
  Vec normal;
  if (intersectionType == L1_WITH_L2) {
    Vec v = Vec_makeFromLine(*l2);
    face = Vec_normalize(v);
  } else {
    Vec v = Vec_makeFromLine(*l1);
    face = Vec_normalize(v);
  }
  normal = Vec_orthogonal(face);

  // Obtain each line's velocity components with respect to the collision
  // face/normal vectors.
  double v1Face = Vec_dotProduct(l1->velocity, face);
  double v2Face = Vec_dotProduct(l2->velocity, face);
  double v1Normal = Vec_dotProduct(l1->velocity, normal);
  double v2Normal = Vec_dotProduct(l2->velocity, normal);

  // Compute the mass of each line (we simply use its length).
  double m1 = Vec_length(Vec_subtract(l1->p1, l1->p2));
  double m2 = Vec_length(Vec_subtract(l2->p1, l2->p2));

  // Perform the collision calculation (computes the new velocities along
  // the direction normal to the collision face such that momentum and
  // kinetic energy are conserved).
  double newV1Normal = ((m1 - m2) / (m1 + m2)) * v1Normal
      + (2 * m2 / (m1 + m2)) * v2Normal;
  double newV2Normal = (2 * m1 / (m1 + m2)) * v1Normal
      + ((m2 - m1) / (m2 + m1)) * v2Normal;

  // Combine the resulting velocities.
  l1->velocity = Vec_add(Vec_multiply(normal, newV1Normal),
                         Vec_multiply(face, v1Face));
  l2->velocity = Vec_add(Vec_multiply(normal, newV2Normal),
                         Vec_multiply(face, v2Face));

  return;
}
