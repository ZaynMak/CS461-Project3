function RayTracer(canvasID, eye, center, up, fov) {
  /**
   * Initializes an instance of the RayTracer class.
   * You may also wish to set up your camera here.
   * (feel free to modify input parameters to set up camera).
   * 
   * @param canvasID (string) - id of the canvas in the DOM where we want to render our image
   */
  // setup the canvas
  this.canvas = document.getElementById(canvasID);

	// setup the camera
	this.camera = new Camera(eye, center, up, fov, this.canvas.width / this.canvas.height);
	

  // setup the background style: current options are 'daylight' or 'white'
  this.sky = 'daylight';

  // initialize the objects and lights
  this.objects = new Array();
  this.lights  = new Array();
	this.triangles = new Array();
}

RayTracer.prototype.draw = function() {
  /**
   * Renders the scene to the canvas.
   * Loops through all pixels and computes a pixel sample at each pixel midpoint.
   * Pixel color should then be computed and assigned to the image.
  **/
  // get the canvas and the image data we will write to
  let context = this.canvas.getContext('2d');
  let image = context.createImageData(this.canvas.width,this.canvas.height);

  // numbers of pixels in x- and y- directions
  const nx = image.width;
  const ny = image.height;
	
  // loop through the canvas pixels
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {

      // compute pixel coordinates in [0,1] x [0,1]
      let px = (i + 0.5) / nx;      // sample at pixel center
      let py = (ny - j - 0.5) / ny; // canvas has y pointing down, but image plane has y going up

			// compute pixel color
      let color = vec3.create();
			
      // YOU NEED TO DETERMINE PIXEL COLOR HERE (see notes)
      // i.e. cast a ray through (px,py) and call some 'color' function
			let r = this.camera.getRayDirection(px, py);
			
			let ray = new Ray(this.camera.eye, vec3.clone(r));

			color = this.color(ray, 2);
			

      // set the pixel color into our final image
      this.setPixel( image , i,j , color[0] , color[1] , color[2] );
    }
  }
  context.putImageData(image,0,0);
}

RayTracer.prototype.color = function(ray, depth) {
	let current_min = new Hit_Info(undefined, undefined, Infinity);
	let curr_index = 0;
	let hit = false;
	
	if (depth == 0){
		return this.background(ray);
	}

	if (this.triangles.length != 0){
		for (let a = 0; a < this.triangles.length / 3; a++){
			let m = intersectTriangle(ray, this.triangles[3*a], this.triangles[3*a]+1, this.triangles[3*a+2], this.camera.eye);
		}
	}

	for (let k = 0; k < this.objects.length; k++){
		let intersection_info = this.objects[k].intersect( ray,  0.001, 1e20);
		if (intersection_info == undefined){
			color = this.background(ray);
		} else {
			hit = true;
			
			// saves in current min the closest intersection point's info
			if (intersection_info.ray_value < current_min.ray_value){
				current_min = intersection_info;
				curr_index = k;
			}
		}
	}

	if (hit){
		// if statement for checkpoint 1 and 2
		if (this.objects[curr_index].material == undefined){
			color = this.objects[curr_index].color;

		// else statement for checkpoint 3 and 4
		} else {
			// calculating la for ambient light
			let la = vec3.fromValues(0.0, 0.0, 0.0);
			for (let l = 0; l < this.lights.length; l++){
				vec3.add(la, la, this.lights[l].color);
			}
			vec3.scale(la, la, 1/this.lights.length);

			// initialize color from the ambient light
			vec3.multiply(color, this.objects[curr_index].material.ka, la);

			// calculate diffuse and specular color
			for (let l = 0; l < this.lights.length; l++){
				let direction_to_light = vec3.subtract(vec3.create(), this.lights[l].location, current_min.intersection_coord);
				vec3.normalize(direction_to_light, direction_to_light);
				
				let shadow_ray = new Ray(current_min.intersection_coord, direction_to_light);
				for (let k = 0; k < this.objects.length; k++){
					let blocking_object = this.objects[k].intersect(shadow_ray, 0.001, 1e20 );

					// if no objects are blocking the light, add shade
					if (blocking_object != undefined){
						continue;
					} else {
						// let shade = this.objects[k].material.shade(ray, this.lights[l], current_min); // k or curr_index
						let shade = this.objects[curr_index].material.shade(ray, this.lights[l], current_min);
						vec3.add(color, color, shade);
					}
				}
			}

			let scattered_ray = new Ray(current_min.intersection_coord, this.objects[curr_index].material.scatter(ray, current_min));

			if (scattered_ray.direction == undefined){
				return color
			}

			let color_scattered = this.color(scattered_ray, depth-1);

			// let scatter_mix = vec3.scale(vec3.create(), color_scattered, 0.3);
			// vec3.scaleAndAdd(color, scatter_mix, color, 0.7);
			

			let scatter_mix = vec3.scale(vec3.create(), color_scattered, this.objects[curr_index].material.weight);
			vec3.scaleAndAdd(color, scatter_mix, color, 1-this.objects[curr_index].material.weight);

			return color;
		}
	}
	return color;
}

RayTracer.prototype.background = function(ray) {
  /**
   * Computes the background color for a ray that goes off into the distance.
   * 
   * @param ray - ray with a 'direction' (vec3) and 'point' (vec3)
   * @returns a color as a vec3 with (r,g,b) values within [0,1]
   * 
   * Note: this assumes a Ray class that has member variable ray.direction.
   * If you change the name of this member variable, then change the ray.direction[1] accordingly.
  **/
  if (this.sky === 'white') {
    // a white sky
    return vec3.fromValues(1,1,1);
  }
  else if (this.sky === 'daylight') {
    // a light blue sky :)
    let t = 0.5*ray.direction[1] + 0.2; // uses the y-values of ray.direction
    if (ray.direction == undefined) t = 0.2; // remove this if you have a different name for ray.direction
    let color = vec3.create();
    vec3.lerp( color , vec3.fromValues(.5,.7,1.)  , vec3.fromValues(1,1,1) , t );
    return color;
  }
  else
    alert('unknown sky ',this.sky);
}

RayTracer.prototype.setPixel = function( image , x , y , r , g , b ) {
  /**
   * Sets the pixel color into the image data that is ultimately shown on the canvas.
   * 
   * @param image - image data to write to
   * @param x,y - pixel coordinates within [0,0] x [canvas.width,canvas.height]
   * @param r,g,b - color to assign to pixel, each channel is within [0,1]
   * @returns none
   * 
   * You do not need to change this function.
  **/
  let offset = (image.width * y + x) * 4;
  image.data[offset  ] = 255*Math.min(r,1.0);
  image.data[offset+1] = 255*Math.min(g,1.0);
  image.data[offset+2] = 255*Math.min(b,1.0);
  image.data[offset+3] = 255; // alpha: transparent [0-255] opaque
}

function Sphere(params) {
  // represents a sphere object
  this.center   = params['center']; // center of the sphere (vec3)
  this.radius   = params['radius']; // radius of sphere (float)
  this.material = params['material']; // material used to shade the sphere (see 'Material' below)
	this.color   = params['color']; // radius of sphere (float)
  this.name     = params['name'] || 'sphere'; // a name to identify the sphere (useful for debugging) (string)
}

Sphere.prototype.intersect = function (ray , tmin , tmax){
	let xo = vec3.create();
	vec3.sub(xo, ray.point, this.center);
  let B = vec3.dot( ray.direction , xo );
  // let C = vec3.dot(xo, xo) - (this.radius * this.radius);
	let C = vec3.sqrDist(ray.point, this.center) - (this.radius * this.radius);
  let discriminant = B*B - C;
  
	if (discriminant < 0.0){
    return undefined;
  } 
	let t = 0;

  let t1 = -B - Math.sqrt(discriminant);
	let t2 = -B + Math.sqrt(discriminant);
	if (t1 < t2) {
		t = t1;
	} else {
		t = t2;
	}
	if ((t < tmax) && (t > tmin)){
		let intersection_coord = vec3.scaleAndAdd(vec3.create(), ray.point, ray.direction, t);
		let normal = vec3.subtract(vec3.create(), intersection_coord, this.center);
		vec3.normalize(normal, normal);
		let hit_info = new Hit_Info(intersection_coord, normal, t);
		return hit_info;
	} else {
		return undefined;
	}
}

function Light(params) {
  // describes a point light source, storing the location of the light
  // as well as ambient, diffuse and specular components of the light
  this.location = params.location; // location of the light
  this.color    = params.color || vec3.fromValues(1,1,1); // default to white (vec3)
  // you might also want to save some La, Ld, Ls and/or compute these from 'this.color'
	this.ld = vec3.clone(this.color);
	this.ls = vec3.clone(this.color);
}

function Material( params ) {
  // represents a generic material class
  this.type  = params.type; // diffuse, reflective, refractive (string)
  this.shine = params.shine; // phong exponent (float)
	this.refractive_index = params.refractive_index; // (float)
	this.weight = params.weight; // float
  this.color = params.color || vec3.fromValues(0.5,0.5,0.5); // default to gray color (vec3)
  // you might also want to save some ka, kd, ks and/or compute these from 'this.color'
	this.ka = vec3.scale(vec3.create(), this.color, 0.5);
	this.kd = vec3.clone(this.color);
	vec3.scale(this.kd, this.kd, 0.3);
	this.ks = vec3.clone(this.color);
}

Material.prototype.shade = function(ray, light, hit_info){
	// calculating diffuse light
	let normal = hit_info.normal;
	let direction_light = vec3.subtract(vec3.create(), light.location, hit_info.intersection_coord);
	vec3.normalize(direction_light, direction_light);

	let kld = vec3.multiply(vec3.create(), this.kd, light.ld);
	let ndotl = vec3.dot(direction_light, normal);
	let max = Math.max(0.0, ndotl);
	vec3.scale(kld, kld, max);

	// if shiny, add the specular light
	let kls = vec3.fromValues(0.0,0.0,0.0);
	if(this.shine != undefined){
		kls = vec3.multiply(vec3.create(), this.ks, light.ls);
		let r =	vec3.scale(vec3.create(), normal, 2 * ndotl);
		vec3.subtract(r, r, direction_light);
		vec3.normalize(r, r);
		let ndotr = vec3.dot(normal, r);
		ndotr = Math.pow(ndotr, this.shine);
		max = Math.max(0.0, ndotr);
		vec3.scale(kls, kls, max);
	}
		
	return vec3.add(vec3.create(), kld, kls);
}

Material.prototype.scatter = function(ray, hit_info){
	let r = undefined;
	if (this.type == 'reflective'){
		r = reflect(ray.direction, hit_info.normal);
	} else if (this.type == 'refractive'){
		r = refract(ray.direction, hit_info.normal, 1/this.refractive_index); // air over refractive index
	}
	return r;
}

//(vec3 eye, vec3 center, vec3 up, float fov, float aspect)
function Camera (eye, center, up, fov, aspect) {
	this.eye = eye;
	this.center = center;
	this.up = up;
	this.fov = fov;
	this.aspect = aspect;
}

Camera.prototype.getRayDirection = function(px, py){
	let d = vec3.distance(this.eye, this.center);

	let scaled_fov = this.fov / 2;

	let h = 2.0 * d * Math.tan(scaled_fov);
	let width = this.aspect * h;
  
	let g = vec3.sub(vec3.create(), this.center, this.eye);
	let normalizedg = vec3.normalize(vec3.create(), g);
	
	let t = vec3.clone(this.up);

	let gcrosst = vec3.cross(vec3.create(), g, t);
	let u = vec3.normalize(vec3.create(), gcrosst);
	let w = vec3.negate(vec3.create(), normalizedg);
  let v = vec3.cross(vec3.create(), w, u);

	let pu = -width/2.0 + width*px;
	let pv = -h/2.0 + h*py;
	let pw = -d;

	let r = vec3.create();
	vec3.scale(r, u, pu);
	vec3.scaleAndAdd(r, r, v, pv);
	vec3.scaleAndAdd(r, r, w, pw);
	let final_r = vec3.create();
	vec3.normalize(final_r, r);
	return final_r
}

function Ray (point, direction){
	this.point = point;
	this.direction = direction;
}

function Hit_Info (intersection_coord, normal, ray_value) {
	this.intersection_coord = intersection_coord;
	this.normal = normal;
	this.ray_value = ray_value;
}

reflect = function( v , n ) {
  /**
    @param v - vec3 for the vector we want to reflect (assume unit vector)
    @param n - vec3 for the vector about which we reflect (unit normal vector)
    @return vec3 for the reflection direction
  **/
  let r = vec3.create();
  vec3.scale( r , n , -2.0*vec3.dot(v,n) );
  vec3.add( r , r , v );
  return r;
}

refract = function( v , n , n1_over_n2 ) {
  /**
   @param n1_over_n2 - (float) ratio of refractive indices of material 1 (where ray comes from) and material 2 (where ray goes into)
   @param v - (vec3) ray direction pointing from material 1 into material 2
   @param n - (vec3) vector normal to the surface which points INTO material 1
  **/
  let dt = vec3.dot(v,n);

	if (dt > 0.0){
		vec3.scale(n,n,-1);
		n1_over_n2 = 1 / n1_over_n2;
	}
	dt = vec3.dot(v,n);
	
  let discriminant = 1.0 - n1_over_n2*n1_over_n2*( 1.0 - dt*dt );
  if (discriminant < 0.0) {
    // total internal reflection
    return reflect(v,n);
  }

  let r1 = vec3.create();
  vec3.scaleAndAdd( r1 , v , n , -dt );
  vec3.scale( r1 , r1 , n1_over_n2 );

  let r2 = vec3.create();
  vec3.scale( r2 , n , Math.sqrt(discriminant) );

  let r = vec3.create();
  vec3.subtract(r,r1,r2);
  return r;
}

let intersectTriangle = function(r, p1, p2, p3, eye) {
  /**
   * takes in a ray direction 'r' and uses the three triangle
   * points defined above to determine the intersection of the 
   * ray with the triangle
   * 
   * @param: r (vec3) ray direction
   * @return: coordinates of intersection point (vec3) or 'undefined' if no intersection
   * 
   * hints: see gl-matrix documentation for 'invert' and 'transformMat3'
  **/
  let A = mat3.create();
  let b = vec3.create();

  for (let d = 0; d < 3; d++) {
    A[d]   = p1[d] - p3[d];
    A[3+d] = p2[d] - p3[d];
    A[6+d] = -r[d];
    b[d]   = eye[d] - p3[d];
  }

  let Ainv = mat3.create();
  mat3.invert( Ainv , A );

  let c = vec3.create();
  vec3.transformMat3( c , b , Ainv );

  let alpha = c[0];
  let beta = c[1];
  let t = c[2];

  let gamma = 1. - alpha - beta;
  if (alpha < 0.0 || alpha > 1.0) return undefined;
  if (beta < 0.0 || beta > 1.0) return undefined;
  if (gamma < 0.0 || gamma > 1.0) return undefined;

  let p = vec3.create();
  for (let d = 0; d < 3; d++)
    p[d] = eye[d] + t * r[d];

  return p;
}