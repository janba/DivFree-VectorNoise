/*

  Accompanying code, Improving Curl Noise, SIGGRAPH Asia 2025
  ---

  This code warps a 2D image using divergence-free vector noise.

*/

//CURL NOISE CONTROLS

//Use anti aliasing (set to one to enable - very slow)
#define AA 0

// Use reprojection (set to zero to disable).
#define REPROJECT 1

// Use 4th order Runge-Kutta (Euler if set to zero)
#define USE_RK4 1

// Number of integration steps. 
// SET TO 2048 FOR PAPER IDENTICAL RESULTS (slow though)
const int N = 64;


// The MIT License
// Copyright © 2017 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// https://www.youtube.com/c/InigoQuilez
// https://iquilezles.org/

// Computes the analytic derivatives of a 3D Gradient Noise. This can be used for example to compute normals to a
// 3d rocks based on Gradient Noise without approximating the gradient by having to take central differences.
//
// More info here: https://iquilezles.org/articles/gradientnoise

// All noise functions here:
//
// https://www.shadertoy.com/playlist/fXlXzf&from=0&num=12


// 0: integer hash
// 1: float hash (aliasing based)
#define METHOD 0


// 0: cubic
// 1: quintic
#define INTERPOLANT 1


#if METHOD==0
vec3 hash( ivec3 p )     // this hash is not production ready, please
{                        // replace this by something better
	ivec3 n = ivec3( p.x*127 + p.y*311 + p.z*74,
                     p.x*269 + p.y*183 + p.z*246,
                     p.x*113 + p.y*271 + p.z*124);

	// 1D hash by Hugo Elias
	n = (n << 13) ^ n;
    n = n * (n * n * 15731 + 789221) + 1376312589;
    return -1.0+2.0*vec3( n & ivec3(0x0fffffff))/float(0x0fffffff);
}
#else
vec3 hash( vec3 p )      // this hash is not production ready, please
{                        // replace this by something better
	p = vec3( dot(p,vec3(127.1,311.7, 74.7)),
			  dot(p,vec3(269.5,183.3,246.1)),
			  dot(p,vec3(113.5,271.9,124.6)));

	return -1.0 + 2.0*fract(sin(p)*43758.5453123);
}
#endif


// return value noise (in x) and its derivatives (in yzw)
vec4 noised( in vec3 x )
{
    // grid
    #if METHOD==0
    ivec3 i = ivec3(floor(x));
    #else
    vec3 i = floor(x);
    #endif
    vec3 f = fract(x);
    
    #if INTERPOLANT==1
    // quintic interpolant
    vec3 u = f*f*f*(f*(f*6.0-15.0)+10.0);
    vec3 du = 30.0*f*f*(f*(f-2.0)+1.0);
    #else
    // cubic interpolant
    vec3 u = f*f*(3.0-2.0*f);
    vec3 du = 6.0*f*(1.0-f);
    #endif    
    
    // gradients
    #if METHOD==0
    vec3 ga = hash( i+ivec3(0,0,0) );
    vec3 gb = hash( i+ivec3(1,0,0) );
    vec3 gc = hash( i+ivec3(0,1,0) );
    vec3 gd = hash( i+ivec3(1,1,0) );
    vec3 ge = hash( i+ivec3(0,0,1) );
	vec3 gf = hash( i+ivec3(1,0,1) );
    vec3 gg = hash( i+ivec3(0,1,1) );
    vec3 gh = hash( i+ivec3(1,1,1) );
    #else
    vec3 ga = hash( i+vec3(0.0,0.0,0.0) );
    vec3 gb = hash( i+vec3(1.0,0.0,0.0) );
    vec3 gc = hash( i+vec3(0.0,1.0,0.0) );
    vec3 gd = hash( i+vec3(1.0,1.0,0.0) );
    vec3 ge = hash( i+vec3(0.0,0.0,1.0) );
	vec3 gf = hash( i+vec3(1.0,0.0,1.0) );
    vec3 gg = hash( i+vec3(0.0,1.0,1.0) );
    vec3 gh = hash( i+vec3(1.0,1.0,1.0) );
    #endif
    
    // projections
    float va = dot( ga, f-vec3(0.0,0.0,0.0) );
    float vb = dot( gb, f-vec3(1.0,0.0,0.0) );
    float vc = dot( gc, f-vec3(0.0,1.0,0.0) );
    float vd = dot( gd, f-vec3(1.0,1.0,0.0) );
    float ve = dot( ge, f-vec3(0.0,0.0,1.0) );
    float vf = dot( gf, f-vec3(1.0,0.0,1.0) );
    float vg = dot( gg, f-vec3(0.0,1.0,1.0) );
    float vh = dot( gh, f-vec3(1.0,1.0,1.0) );
	
    // interpolations
    return vec4( va + u.x*(vb-va) + u.y*(vc-va) + u.z*(ve-va) + u.x*u.y*(va-vb-vc+vd) + u.y*u.z*(va-vc-ve+vg) + u.z*u.x*(va-vb-ve+vf) + (-va+vb+vc-vd+ve-vf-vg+vh)*u.x*u.y*u.z,    // value
                 ga + u.x*(gb-ga) + u.y*(gc-ga) + u.z*(ge-ga) + u.x*u.y*(ga-gb-gc+gd) + u.y*u.z*(ga-gc-ge+gg) + u.z*u.x*(ga-gb-ge+gf) + (-ga+gb+gc-gd+ge-gf-gg+gh)*u.x*u.y*u.z +   // derivatives
                 du * (vec3(vb,vc,ve) - va + u.yzx*vec3(va-vb-vc+vd,va-vc-ve+vg,va-vb-ve+vf) + u.zxy*vec3(va-vb-ve+vf,va-vb-vc+vd,va-vc-ve+vg) + u.yzx*u.zxy*(-va+vb+vc-vd+ve-vf-vg+vh) ));
}


// Number of frequencies to use for Curl Noise. Increasing 
// this number can improve quality, but it also decreases fps
float SCALE = 1.0;

const float PI = 3.141592653589793;

// Step length used for streamline tracing (t)
// Increase this to make point placement more random, decrease to make regular.
float step_len = 1.1; // Useful range approx [0-4]

float noise_scale = 7.0;


// Returns a vector in the Curl Noise vector field
vec3 G(in vec3 x) {
    return noised(x*noise_scale).yzw*noise_scale;
}

// Returns a vector in the Curl Noise vector field
vec3 V(in vec3 x, in vec3 Z) {
    vec3 g = noised(x*noise_scale).yzw*noise_scale;
    return cross(g, Z);
}


float f(in vec3 x) {
    vec4 g = noised(x*noise_scale);
    return g.x;
}

const vec3 Z = vec3(0,0,1);

vec2 simple(in vec2 xp, in mat3 R, float t) {
    vec3 x = vec3(xp, 0.0);
    float tau = noised(x*noise_scale).x;
    float dt = t / float(N);
    for (int i=0;i<N;++i) {
        vec3 v = cross(noised(x*noise_scale).yzw*noise_scale, Z);
        x = x + dt * v;
#if REPROJECT != 0
        for (int j=0; j<REPROJECT; ++j) {
            vec4 n = noised(x*noise_scale);
            vec3 g = n.yzw*noise_scale;
            x = x - g*(n.x-tau)/(dot(g,g)+1e-30);
            x -= dot(x,Z) * Z;
        }
#endif
    }
    return x.xy;
}


// Use fourth order Runge-Kutta to trace a streamline for a given
// time step. Note that since the vector field is the curl noise
// vector field, this function is where the actual jittering takes
// place.
vec2 RK4(in vec2 xp, float t) {
    vec3 x = vec3(xp, 0.0);
    float tau = f(x);
    float dt = t / float(N);
    for (int i=0;i<N;++i) {
        vec3 a = dt * V(x,Z);
        vec3 b = dt * V(x+a/2.0,Z);
        vec3 c = dt * V(x+b/2.0,Z);
        vec3 d = dt * V(x+c,Z);
        x = x + (a+2.0*b+2.0*c+d)/6.0;
#if REPROJECT!=0
        for (int j=0; j<REPROJECT; ++j) {
            vec4 n = noised(x*noise_scale);
            vec3 g = n.yzw*noise_scale;
            x = x - g*(n.x-tau)/(dot(g,g)+1e-30);
            x -= dot(x,Z) * Z;
        }
#endif
    }
    return x.xy;
}



vec2 uvmap(vec2 uv) {
    return vec2(abs(uv.x), uv.y + 0.05*cos(uv.x*12.0*PI));
}


float saturate( float x ) { return clamp( x, 0.0, 1.0 ); }

// The color mapping
vec3 viridis_quintic( float x )
{
	x = saturate( x );
	vec4 x1 = vec4( 1.0, x, x * x, x * x * x ); // 1 x x2 x3
	vec4 x2 = x1 * x1.w * x; // x4 x5 x6 x7
	return vec3(
		dot( x1.xyzw, vec4( +0.280268003, -0.143510503, +2.225793877, -14.815088879 ) ) + dot( x2.xy, vec2( +25.212752309, -11.772589584 ) ),
		dot( x1.xyzw, vec4( -0.002117546, +1.617109353, -1.909305070, +2.701152864 ) ) + dot( x2.xy, vec2( -1.685288385, +0.178738871 ) ),
		dot( x1.xyzw, vec4( +0.300805501, +2.614650302, -12.019139090, +28.933559110 ) ) + dot( x2.xy, vec2( -33.491294770, +13.762053843 ) ) );
}


mat3 rot_mat_Y(float t) {
    float cosT = cos(t);
    float sinT = sin(t);
    return mat3(cosT, 0.0, sinT, 
                0.0, 1.0, 0.0,-
                sinT, 0.0, cosT);
}

float f_hash(float x) {
    return float((17*int(10.0*mod(x, 1.0)))%13)/13.0;
}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
#if AA==1
    float iTime_s = iTime;// * 0.025;
    for (int x=0;x<4;++x) for (int y=0;y<4;++y) {
    vec2 fragCoord_s = fragCoord;
    fragCoord_s.x += float(x)/4.0;
    fragCoord_s.y += float(y)/4.0;
    
    float t = 0.5*sin(iTime_s * 0.03);// 0.5 USED FOR EXPERIMENTS
    float p = float(int(iTime_s*0.25/(2.0*PI)));
    float scale = SCALE;
    vec2 uv_o = scale*(fragCoord_s/iResolution.x-0.5);
#else
    float t = 0.5*sin(iTime * 0.03);// 0.5 USED FOR EXPERIMENTS
    float p = float(int(iTime*0.25/(2.0*PI)));
    float scale = SCALE;
    vec2 uv_o = scale*(fragCoord/iResolution.x-0.5);
#endif
    vec2 uv  = uv_o;
    mat3 R = rot_mat_Y(t);
#if USE_RK4==1
        uv = RK4(uv, t);
#else
        uv = simple(uv, R, t);
#endif
    uv.y += 0.05*abs(cos(1.0*PI*10.0*(uv.x-uv.y)));
    float n = dot(uv, vec2(0,1))+ 1213.12;
#if AA==1
    fragColor.rgb += 0.0625*viridis_quintic(f_hash(n));
    }
#else
    fragColor.rgb += viridis_quintic(f_hash(n));
#endif

}