#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <arm_neon.h>

#define N 1000
#define STEPS 16

extern void enable_runfast();

float vx[N], vy[N], vz[N], xnew[N], ynew[N], znew[N];
float32_t m[N], x[N], y[N], z[N];

void  diff(struct timespec * difference, struct timespec start, struct timespec end)
{
  if ((end.tv_nsec-start.tv_nsec)<0) {
    difference->tv_sec = end.tv_sec-start.tv_sec-1;
    difference->tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
//    printf("add big number\n");
  } else {
    difference->tv_sec = end.tv_sec-start.tv_sec;
    difference->tv_nsec = end.tv_nsec-start.tv_nsec;
//    printf("no add\n");
//    printf("start: %d s %d ns\n", (int)start.tv_sec, (int)start.tv_nsec);
//    printf("end:   %d s %d ns\n", (int)end.tv_sec, (int)end.tv_nsec);
  }
}

void init(void) {
  int i;

  for(i=0; i<N; i++) { /* Foreach particle "i" ... */
    x[i] = rand();
    y[i] = rand();
    z[i] = rand();
    vx[i] = rand()/100;
    vy[i] = rand()/100;
    vz[i] = rand()/100;
    m[i] = rand();
  }
}

int main (int argc, char * argv[]) {
  int s,i,j;
  float dt=0.001;
  float eps=0.0000001;
  struct timespec t1, t2, d;
  FILE *fp;
  char *outputFilename = "results.txt";
  
  float32x4_t vec_dx, vec_dy, vec_dz;
  float32x4_t vec_dxi, vec_dyi, vec_dzi;
  float32x4_t vec_invr, vec_invr3, vec_f;
//  float32x4_t vec_eps;
  float32x4_t vec_ax, vec_ay, vec_az;
  float32_t ax, ay, az;
  
  //eps
  //vec_eps = vld4_dup_f32(eps);
/*float32x4_t vdupq_n_f32 (float32_t)
  Form of expected instruction(s): vdup.32 q0, r0*/

  enable_runfast();
  init();

  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
//printf("time: %d s %d ns\n", (int)t1.tv_sec, (int)t1.tv_nsec);

  for (s=0; s<STEPS; s++) {
    for(i=0; i<N; i+=4) { /* Foreach particle "i" ... */
      ax=0.0f;
      ay=0.0f;
      az=0.0f;
      vec_dxi = vld1q_f32(&x[i]);
      vec_dyi = vld1q_f32(&y[i]);
      vec_dzi = vld1q_f32(&z[i]);
/*float32x4_t vdupq_n_f32 (float32_t)
  Form of expected instruction(s): vdup.32 q0, r0*/
      for(j=0; j<N; j+=4) { /* Loop over all particles "j" */
	      //dx=x[j]-x[i];
	      //dy=y[j]-y[i];
	      //dz=z[j]-z[i];
	      vec_dx = vld1q_f32(&x[j]);
	      vec_dy = vld1q_f32(&y[j]);
	      vec_dz = vld1q_f32(&z[j]);
	      vec_dx = vsubq_f32(vec_dx, vec_dxi);
	      vec_dy = vsubq_f32(vec_dy, vec_dyi);
	      vec_dz = vsubq_f32(vec_dz, vec_dzi);
/*float32x4_t vld1q_f32 (const float32_t *)
  Form of expected instruction(s): vld1.32 {d0, d1}, [r0]*/
/*float32x4_t vsubq_f32 (float32x4_t, float32x4_t)
  Form of expected instruction(s): vsub.f32 q0, q0, q0*/
	      
	      //invr = 1.0f/sqrtf(dx*dx + dy*dy + dz*dz + eps);
	      vec_invr = vdupq_n_f32(eps);
	      vec_invr = vmlaq_f32(vec_invr, vec_dx, vec_dx);
	      vec_invr = vmlaq_f32(vec_invr, vec_dy, vec_dy);
	      vec_invr = vmlaq_f32(vec_invr, vec_dz, vec_dz);
	      vec_invr = vrsqrteq_f32(vec_invr);
/*float32x4_t vrsqrtsq_f32 (float32x4_t, float32x4_t)
  Form of expected instruction(s): vrsqrts.f32 q0, q0, q0*/
/*float32x4_t vmlaq_f32 (float32x4_t, float32x4_t, float32x4_t)
  Form of expected instruction(s): vmla.f32 q0, q0, q0*/

	      //invr3 = invr*invr*invr;
	      vec_invr3 = vmulq_f32(vec_invr, vec_invr);
	      vec_invr3 = vmulq_f32(vec_invr3, vec_invr);
	      //f=m[j]*invr3;
	      vec_f = vld1q_f32(&m[j]);
	      vec_f = vmulq_f32(vec_f, vec_invr3);
/*float32x4_t vmulq_f32 (float32x4_t, float32x4_t)
  Form of expected instruction(s): vmul.f32 q0, q0, q0*/
        
	      //ax += f*dx; /* accumulate the acceleration from gravitational attraction */
	      //ay += f*dy;
	      //az += f*dx;
	      vec_ax = vmulq_f32(vec_f, vec_dx);
	      vec_ay = vmulq_f32(vec_f, vec_dy);
	      vec_az = vmulq_f32(vec_f, vec_dx);
	      ax += vgetq_lane_f32(vec_ax, 0);
	      ax += vgetq_lane_f32(vec_ax, 1);
	      ax += vgetq_lane_f32(vec_ax, 2);
	      ax += vgetq_lane_f32(vec_ax, 3);
	      
	      ay += vgetq_lane_f32(vec_ay, 0);
	      ay += vgetq_lane_f32(vec_ay, 1);
	      ay += vgetq_lane_f32(vec_ay, 2);
	      ay += vgetq_lane_f32(vec_ay, 3);
	      
	      az += vgetq_lane_f32(vec_az, 0);
	      az += vgetq_lane_f32(vec_az, 1);
	      az += vgetq_lane_f32(vec_az, 2);
	      az += vgetq_lane_f32(vec_az, 3);
	      
/*float32x4_t vmlaq_f32 (float32x4_t, float32x4_t, float32x4_t)
  Form of expected instruction(s): vmla.f32 q0, q0, q0*/
      }
      xnew[i] = x[i] + dt*vx[i] + 0.5f*dt*dt*ax; /* update position of particle "i" */
      ynew[i] = y[i] + dt*vy[i] + 0.5f*dt*dt*ay;
      znew[i] = z[i] + dt*vz[i] + 0.5f*dt*dt*az;
      vx[i] += dt*ax; /* update velocity of particle "i" */
      vy[i] += dt*ay;
      vz[i] += dt*az;
      xnew[i+1] = x[i+1] + dt*vx[i+1] + 0.5f*dt*dt*ax; /* update position of particle "i" */
      ynew[i+1] = y[i+1] + dt*vy[i+1] + 0.5f*dt*dt*ay;
      znew[i+1] = z[i+1] + dt*vz[i+1] + 0.5f*dt*dt*az;
      vx[i+1] += dt*ax; /* update velocity of particle "i" */
      vy[i+1] += dt*ay;
      vz[i+1] += dt*az;
      xnew[i+2] = x[i+2] + dt*vx[i+2] + 0.5f*dt*dt*ax; /* update position of particle "i" */
      ynew[i+2] = y[i+2] + dt*vy[i+2] + 0.5f*dt*dt*ay;
      znew[i+2] = z[i+2] + dt*vz[i+2] + 0.5f*dt*dt*az;
      vx[i+2] += dt*ax; /* update velocity of particle "i" */
      vy[i+2] += dt*ay;
      vz[i+2] += dt*az;
      xnew[i+3] = x[i+3] + dt*vx[i+3] + 0.5f*dt*dt*ax; /* update position of particle "i" */
      ynew[i+3] = y[i+3] + dt*vy[i+3] + 0.5f*dt*dt*ay;
      znew[i+3] = z[i+3] + dt*vz[i+3] + 0.5f*dt*dt*az;
      vx[i+3] += dt*ax; /* update velocity of particle "i" */
      vy[i+3] += dt*ay;
      vz[i+3] += dt*az;
    }
    for(i=0;i<N;i++) { /* copy updated positions back into original arrays */
      x[i] = xnew[i];
      y[i] = ynew[i];
      z[i] = znew[i];
    }
  }
//printf("time: %d s %d ns\n", (int)t1.tv_sec, (int)t1.tv_nsec);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2);
  
  // Print results to file so we can compare to ensure optimizations do not alter functionality
  fp = fopen(outputFilename, "w");
  if (fp == NULL) {
    fprintf(stderr, "Can't open output file %s!\n", outputFilename);
    exit(1);
  }
  for (i=0; i<N; i++) {
	fprintf(fp, "%f %f %f\n", x[i], y[i], z[i]);
  }
  fclose(fp);
  
  diff(&d, t1, t2);
  printf("Execution Time: %ld sec, %ld nsec\n", d.tv_sec, d.tv_nsec);
  return 0;
}
