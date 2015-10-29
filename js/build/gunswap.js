(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
// Most of the ideas below are from http://www.benknowscode.com/2012/09/path-interpolation-using-cubic-bezier_9742.html
// The only thing I can take any credit for is modifying the algorithm to accomodate variable input control points and
// the "matchVelocity" idea.

module.exports.interpolateBezierSpline = function(P,t,v_0,v_T,v_0scale,v_Tscale,matchVelocity) {
		
	// t goes from 0 to 1

	/* dwellPath object looks like this */
	/*
		[
			[{x1,y1,z1},{x2,y2,z2},...]
			,[{x1,y1,z1},{x2,y2,z2},...]
		]			
	*/

	var dwellPosition;

	if (t == 0) {
		dwellPosition = P[0];
	} else if (t == 1) {
		dwellPosition = P.last();
	} else {

		/* if P is just one point, duplicate it */
		if (P.length == 1) {
			P.push(P[0]);
		}			

		var C = [{x: P[0].x+v_0.dx*v_0scale, y: P[0].y+v_0.dy*v_0scale, z: P[0].z+v_0.dz*v_0scale},{x: P.last().x-v_T.dx*v_Tscale, y: P.last().y-v_T.dy*v_Tscale, z: P.last().z-v_T.dz*v_Tscale}];
		var eps = .00001;

		var c = [];
		var path = [];

		for (var i = 0; i < P.length-1; i++) {
			var p0 = P[i];
			var p1 = P[i+1];

			var c0, c1;

			if (i == 0) {
				c0 = C[0];
			} else {
				var c1prev = c[c.length-1];
				var c0 = { x: p0.x + (p0.x - c1prev.x), y: p0.y + (p0.y - c1prev.y), z: p0.z + (p0.z - c1prev.z) };
				c.push(c0);
			}

			if (i+1 == P.length-1) {
				c1 = C[1];
			} else {
				var m0 = { x: (P[i].x+P[i+1].x)/2, y: (P[i].y+P[i+1].y)/2, z: (P[i].z+P[i+1].z)/2 };
				var m1 = { x: (P[i+1].x+P[i+2].x)/2, y: (P[i+1].y+P[i+2].y)/2, z: (P[i+1].z+P[i+2].z)/2 };
				var l0 = Math.sqrt( Math.pow(P[i].x - P[i+1].x,2) + Math.pow(P[i].y - P[i+1].y,2) + Math.pow(P[i].z - P[i+1].z,2) );
				var l1 = Math.sqrt( Math.pow(P[i+1].x - P[i+2].x,2) + Math.pow(P[i+1].y - P[i+2].y,2) + Math.pow(P[i+1].z - P[i+2].z,2) );
				var _t = l0/(l0+l1);
				var q = { x: (1-_t)*m0.x + _t*m1.x, y: (1-_t)*m0.y + _t*m1.y, z: (1-_t)*m0.z + _t*m1.z };
				c1 = { x: p1.x + (m0.x-q.x), y: p1.y + (m0.y-q.y), z: p1.z + (m0.z-q.z) };
				c.push(c1);
			}

			for (var _t = 0; _t <= 1+eps; _t += .02) {
				path.push(
					{
						x: Math.pow(1-_t,3)*p0.x + 3*Math.pow(1-_t,2)*_t*c0.x + 3*(1-_t)*Math.pow(_t,2)*c1.x + Math.pow(_t,3)*p1.x,
						y: Math.pow(1-_t,3)*p0.y + 3*Math.pow(1-_t,2)*_t*c0.y + 3*(1-_t)*Math.pow(_t,2)*c1.y + Math.pow(_t,3)*p1.y,
						z: Math.pow(1-_t,3)*p0.z + 3*Math.pow(1-_t,2)*_t*c0.z + 3*(1-_t)*Math.pow(_t,2)*c1.z + Math.pow(_t,3)*p1.z
					}
				);
				if (path.length == 1) {
					path.last().dist = 0;
				} else {
					path.last().dist = path[path.length-2].dist + Math.sqrt(Math.pow(path.last().x-path[path.length-2].x,2)+Math.pow(path.last().y-path[path.length-2].y,2)+Math.pow(path.last().z-path[path.length-2].z,2));
				}
			}
		}

		var dwellPositionIx;
		if (siteswap.matchVelocity) {
			var v_0mag = Math.sqrt(Math.pow(v_0.dx,2)+Math.pow(v_0.dy,2)+Math.pow(v_0.dz,2));
			var v_Tmag = Math.sqrt(Math.pow(v_T.dx,2)+Math.pow(v_T.dy,2)+Math.pow(v_T.dz,2));
			var L = path.last().dist;
			var T = siteswap.dwellDuration;
			var j = (6*T*(v_Tmag + v_0mag) - 12*L)/Math.pow(T,3);
			var a_0 = (v_Tmag - v_0mag - .5*j*Math.pow(T,2))/T;
			var dt = t*T;
			var p_dt = v_0mag*dt + .5*a_0*Math.pow(dt,2) + (1/6)*j*Math.pow(dt,3);

			dwellPositionIx = 0;
			for (var i = 0; i < path.length; i++) {
				if (path[i].dist <= p_dt) {
					dwellPositionIx = i;
				}
			}
		} else {
			dwellPositionIx = Math.floor(t*path.length);
		}
		
		dwellPosition = path[dwellPositionIx < 0 ? 0 : dwellPositionIx];

	}

	return dwellPosition;

}
},{}],2:[function(require,module,exports){
var util = require('./util.js');

function BounceGA(gaConfig, fitnessConfig) {
	/* TO DO - fill in unprovided inputs */
	this.gaConfig = gaConfig;
	this.fitnessConfig = fitnessConfig;

	/* get axes for surfaces, surfaces will always have bottom edge perpendicular to y-axis */
	if (fitnessConfig.axis1 === undefined) {
		for (var i = 0; i < this.fitnessConfig.surfaces.length; i++) {	
			var surface = this.fitnessConfig.surfaces[i];
			util.normalize(surface.normal);
			var axis1;
			if (surface.normal.x == 0 && surface.normal.z == 0) {
				axis1 = {x:1,y:0,z:0};
			} else {
				axis1 = {x:-surface.normal.z,y:0,z:surface.normal.x};
			}
			var axis2 = util.cross(surface.normal,axis1);
			util.normalize(axis1);
			util.normalize(axis2);
			util.multiply(axis1,surface.scale);
			util.multiply(axis2,surface.scale);
			surface.axis1 = axis1;
			surface.axis2 = axis2;
		}
	}

	this.generations = 0;
	this.rankedPopulation = [];	// the current, ranked, population
	this.fittestMembers = []; // array of the fittest members by generation
	this.ableToFindSolution = undefined;

}

BounceGA.prototype.getBouncePath = function(v) {

	var pt = cloneObject(this.fitnessConfig.p0); 
	var vt = cloneObject(v);

	var p = [];
	var velocities = [];
	p.push(cloneObject(pt));
	velocities.push(cloneObject(vt));

	// reset surfaces 
	for (var i = 0; i < this.fitnessConfig.surfaces.length; i++) {
		this.fitnessConfig.surfaces[i].bounces = 0;
		this.fitnessConfig.surfaces[i].colliding = false;
	}

	var totalBounces = 0;
	var actualBounceOrder = [];

	var pTchoices = [];

	for (var t = 0; t <= this.fitnessConfig.maxT; t += this.fitnessConfig.dt) {

		// update position / velocity
		util.add(pt, util.multiply(cloneObject(vt), this.fitnessConfig.dt));
		vt.y += this.fitnessConfig.G*this.fitnessConfig.dt;

		p.push(cloneObject(pt));
		velocities.push(cloneObject(vt));

		// check for collisions
		for (var i = 0; i < this.fitnessConfig.surfaces.length; i++) {
			var surface = this.fitnessConfig.surfaces[i];	
			var distance = (surface.normal.x*pt.x + surface.normal.y*pt.y + surface.normal.z*pt.z - surface.normal.x*surface.position.x - surface.normal.y*surface.position.y - surface.normal.z*surface.position.z);

			if (
				Math.abs(distance) <= this.fitnessConfig.R // distance from plane is within radius
				&& Math.abs(util.dot(util.sub(cloneObject(pt),surface.position),util.normalize(cloneObject(surface.axis1)))) <= surface.scale
				&& Math.abs(util.dot(util.sub(cloneObject(pt),surface.position),util.normalize(cloneObject(surface.axis2)))) <= surface.scale
				&& !surface.colliding
			) {
				surface.colliding = true;
				surface.bounces++;
				actualBounceOrder.push(i);
				totalBounces++;

				//new velocity from bounce
				var bounceV = (new THREE.Vector3()).copy(surface.normal);
				var bounceV = cloneObject(surface.normal);
				util.multiply(bounceV, util.dot(vt, surface.normal));
				util.multiply(bounceV, 2*this.fitnessConfig.C);
				util.negate(bounceV);
				util.add(vt, bounceV);

			} else if (Math.abs(distance) > this.fitnessConfig.R && surface.colliding) {
				surface.colliding = false;
			}

		}

		if (t >= this.fitnessConfig.minT || t + this.fitnessConfig.dt > this.fitnessConfig.maxT) {	
			/* if bounce order is specified, check that it is followed */
			var correctBounceOrder = true;
			if (this.fitnessConfig.surfaceBounceOrder && actualBounceOrder.toString() != this.fitnessConfig.surfaceBounceOrder.toString()) {
				correctBounceOrder = false;
			}

			var enoughBounces = false;
			/* check that the total number of bounces was correct */
			if (totalBounces == this.fitnessConfig.numBounces) {
				enoughBounces = true;
			}

			// must match expected number of bounces, otherwise fitness is terrible
			if ( 
				correctBounceOrder
				&& enoughBounces
				&& ( (this.fitnessConfig.catchSign == 1 && vt.y >= 0) || (this.fitnessConfig.catchSign == -1 && vt.y <= 0 ) || this.fitnessConfig.catchSign == undefined ) 
				&& ( (this.fitnessConfig.tossSign == 1 && v.y >= 0) || (this.fitnessConfig.tossSign == -1 && v.y <= 0 ) || this.fitnessConfig.tossSign == undefined ) 
			) {
				pTchoices.push({
					pT: cloneObject(pt), 
					T: t,
					actualpT: this.fitnessConfig.pT // need this for comparator function below
				});
				//return {path: p, T: t};
			}
		}		

	}

	if (pTchoices.length > 0) {
		// go through pT choices and grab the best
		pTchoices.sort(function(a,b) { 
			return util.magnitude(util.sub(cloneObject(a.pT),a.actualpT)) - util.magnitude(util.sub(cloneObject(b.pT),a.actualpT));
		});

		return {path: p.splice(0,pTchoices[0].T/this.fitnessConfig.dt), velocities: velocities.splice(0,pTchoices[0].T/this.fitnessConfig.dt), T: pTchoices[0].T};
	} else {
		return {path: undefined, velocities: undefined, T: undefined};
	}

}

BounceGA.prototype.generateMember = function(parent) {	

	var v = {x:0, y:0, z:0};

	if (parent && parent.fitness < 100) {

		if (Math.random() < this.gaConfig.mutationChance) {			
			v.x = (1-2*Math.random())*(parent.fitness < this.gaConfig.mutationScale ? parent.fitness : this.gaConfig.mutationScale);
			v.y = (1-2*Math.random())*(parent.fitness < this.gaConfig.mutationScale ? parent.fitness : this.gaConfig.mutationScale);
			v.z = (1-2*Math.random())*(parent.fitness < this.gaConfig.mutationScale ? parent.fitness : this.gaConfig.mutationScale);
			util.add(v, parent.v);
		} else {
			v = cloneObject(parent.v);
		}
	} else {

		// v.x = (1-2*Math.random())*this.gaConfig.initialScale;
		// v.y = (1-2*Math.random())*this.gaConfig.initialScale;
		// v.z = (1-2*Math.random())*this.gaConfig.initialScale;
		
		// should always be tossing towards the first bounce surface
		// TO DO - get this working...

		var targetSurface;
		if(this.fitnessConfig.surfaceBounceOrder) {
			targetSurface = this.fitnessConfig.surfaces[this.fitnessConfig.surfaceBounceOrder[0]];
		} else {
			targetSurface = this.fitnessConfig.surfaces[Math.floor(Math.random()*this.fitnessConfig.surfaces.length)];
		}
		
		// first find a random spot on the surface
		v = cloneObject(targetSurface.position);
		util.add(v, util.multiply(cloneObject(targetSurface.axis1), 1-2*Math.random()));
		util.add(v, util.multiply(cloneObject(targetSurface.axis2), 1-2*Math.random()));

		// save and remove the y component
		var targetY = v.y;
		v.y = 0;

		var p0copy = cloneObject(this.fitnessConfig.p0);
		p0copy.y = 0;
		var pTcopy = cloneObject(this.fitnessConfig.pT);
		pTcopy.y = 0;

		// the minimum possible velocity is for us to get to this spot at the half-way point
		// var minV = p0copy.distanceTo(v)+pTcopy.distanceTo(v) / (this.fitnessConfig.maxT);
		// maximum velocity is for us to go directly to the spot and back in minT
		// ACTUALLY THIS SHOULD BE HIGHER
		// var maxV = (p0copy.distanceTo(v)+pTcopy.distanceTo(v) / (this.fitnessConfig.minT))/this.fitnessConfig.C;
		var minV = 0;
		var maxV = this.gaConfig.initialScale;

		// create velocity vector from p0 to random spot on surface
		util.sub(v, p0copy);
		
		util.normalize(v);
		util.multiply(v, minV + (maxV-minV)*Math.random());

		// now get minV and maxV for y-component
		// the min velocity (or max velocity down) is either 0 if tossSign is 1 or targetY is above, or the time it takes to bounce back from straight down
		// this calculation is naively simplified to not consider C or G, may need to refactor
		minV = (this.fitnessConfig.tossSign == 1 || targetY > this.fitnessConfig.p0.y) ? 0 : -( (this.fitnessConfig.p0.y - targetY)*2 / (this.fitnessConfig.minT) );
		// the max velocity is either 0 if tossSign is -1 or the time it takes to get back from a straight toss up
		//maxV = (this.fitnessConfig.pT.y - this.fitnessConfig.p0.y - .5*this.fitnessConfig.G*this.fitnessConfig.maxT*this.fitnessConfig.maxT)/this.fitnessConfig.maxT;
		maxV = this.gaConfig.initialScale;

		v.y = minV + (maxV-minV)*Math.random();
	
	}

	var bouncePath = this.getBouncePath(v);

	return {
		v: v, 
		T: bouncePath.T,
		fitness: bouncePath.path ? util.magnitude(util.sub(cloneObject(bouncePath.path[bouncePath.path.length-1]), this.fitnessConfig.pT)) : 100
	};
}

BounceGA.prototype.evolve = function() {
	if (this.generations >= this.gaConfig.maxGenerations) {
		this.ableToFindSolution = false;
	} else if (this.rankedPopulation.length > 0 && this.rankedPopulation[0].fitness <= this.gaConfig.fitnessThreshold) {
		this.ableToFindSolution = true;
	} else {
		
		/* if this is first generation, we have no viable solutions yet, or we the noGA flag is set create them all from scratch */
		if (this.rankedPopulation.length == 0 || this.gaConfig.noGA) {
			this.rankedPopulation = []; // just make sure the population is empty
			for (var i = 0; i < this.gaConfig.populationSize; i++) {
				this.rankedPopulation.push(this.generateMember(undefined));
			}				
		} 
		/* otherwise create the next generation */
		else {
			/* we want to create a selection bias based on the fitness */
			var totalWeights = 0;
			var selectorHelper = []; 
			for (var i = 0; i < this.rankedPopulation.length; i++) {
				totalWeights += 1/this.rankedPopulation[i].fitness;
				selectorHelper.push(totalWeights);
			}
			/* select and mutate to create new members using the weighted selector helper array, so we'll pick more of the better members */
			var newPopulation = [];
			for (var i = 0; i < this.gaConfig.populationSize; i++) {
				// make half the population completely new members, the other half choose based on fitness
				if (i < this.gaConfig.populationSize*.5) {
					newPopulation.push(this.generateMember(undefined));
				} else {
					var rand = Math.random()*totalWeights;
					for (var j = 0; j < selectorHelper.length; j++) {
						if (rand < selectorHelper[j]) {
							newPopulation.push(this.generateMember(this.rankedPopulation[j]));
							break;
						}
					}
				}		
			}
			this.rankedPopulation = newPopulation;			
		}

		/* sort population, hence the variable name */
		this.rankedPopulation.sort(function(a,b) { return a.fitness - b.fitness; });
		var fittestMember = this.rankedPopulation[0];
		this.fittestMembers.push(fittestMember);
		this.generations++;

		if (this.generations % 100 == 0) {
			console.log(this.generations + ' ' + fittestMember.fitness);
		}

		/* call evolve again, make sure to use setTimeout so we don't lock the browser */		
		this.evolve();
		/*
		var self = this;
		setTimeout(function() { self.evolve(); }, 0);
		*/
	}
}

module.exports.BounceGA = BounceGA;

// function Animator(ga) {

// 	var 
// 		container,
// 		width,
// 		height,
// 		camera, 
// 		scene, 
// 		renderer,
// 		/* camera starting point */
// 		camTheta = Math.PI, 
// 		camPhi = .3, 
// 		camRadius = 5,
// 		/* helpers for mouse interaction */
// 		isMouseDown = false, 
// 		onMouseDownTheta, 
// 		onMouseDownPhi, 
// 		onMouseDownPosition;

// 	container = $('#animationContainer');

// 	width = container.width();
// 	height = container.height();

// 	camera = new THREE.PerspectiveCamera( 75, width / height, .05, 100 );
// 	updateCamera();

// 	scene = new THREE.Scene();

// 	/* lights */
// 	var ceilingLight = new THREE.PointLight( 0xffffff );
// 	ceilingLight.position.set(0,20,0);
// 	scene.util.add( ceilingLight );

// 	/* surfaces */
// 	ga.fitnessConfig.surfaces.map(function(a) {
// 		var surface = {
// 			position: new THREE.Vector3(a.position.x,a.position.y,a.position.z),
// 			normal: new THREE.Vector3(a.normal.x,a.normal.y,a.normal.z),
// 			axis1: new THREE.Vector3(a.axis1.x,a.axis1.y,a.axis1.z),
// 			axis2: new THREE.Vector3(a.axis2.x,a.axis2.y,a.axis2.z),
// 			scale: a.scale
// 		}
// 		var surfaceGeom = new THREE.Geometry();
// 		surfaceGeom.vertices.push( (new THREE.Vector3()).copy(surface.position).util.add((new THREE.Vector3()).util.add(surface.axis1).util.add(surface.axis2)) );
// 		surfaceGeom.vertices.push( (new THREE.Vector3()).copy(surface.position).util.add((new THREE.Vector3()).util.add(surface.axis1).util.negate().util.add(surface.axis2)) );
// 		surfaceGeom.vertices.push( (new THREE.Vector3()).copy(surface.position).util.add((new THREE.Vector3()).util.add(surface.axis1).util.add(surface.axis2).util.negate()) );
// 		surfaceGeom.vertices.push( (new THREE.Vector3()).copy(surface.position).util.add((new THREE.Vector3()).util.add(surface.axis2).util.negate().util.add(surface.axis1)) );

// 		surfaceGeom.faces.push( new THREE.Face3( 0, 1, 2 ) );
// 		surfaceGeom.faces.push( new THREE.Face3( 2, 0, 3 ) );

// 		var surfaceMesh = new THREE.Mesh(surfaceGeom, new THREE.MeshBasicMaterial( { color: 'grey', side: THREE.DoubleSide } ));
// 		scene.util.add(surfaceMesh);
// 	});

// 	/* draw ball */

// 	var ballMesh = new THREE.Mesh(new THREE.SphereGeometry( .1, 32, 32 ), new THREE.MeshBasicMaterial( { color: 'red' } ));
// 	ballMesh.position = ga.fitnessConfig.p0;
// 	scene.util.add(ballMesh);

// 	var targetMesh = new THREE.Mesh(new THREE.SphereGeometry( .1, 32, 32 ), new THREE.MeshBasicMaterial( { color: 'green' } ));
// 	targetMesh.position = ga.fitnessConfig.pT;
// 	scene.util.add(targetMesh);

// 	/* add axes for debugging */
// 	axes = buildAxes( 1 );
// 	scene.util.add(axes);

// 	/* create the renderer and add it to the canvas container */
// 	if( !window.WebGLRenderingContext ) {
// 		renderer = new THREE.CanvasRenderer();	
// 	} else {
// 		renderer = new THREE.WebGLRenderer( {antialias: true} );
// 	}

// 	renderer.setSize( width, height );

// 	container.empty();
// 	container.append(renderer.domElement);

// 	//add the event listeners for mouse interaction
// 	renderer.domElement.addEventListener( 'mousemove', onDocumentMouseMove, false );
// 	renderer.domElement.addEventListener( 'mousedown', onDocumentMouseDown, false );
// 	renderer.domElement.addEventListener( 'mouseup', onDocumentMouseUp, false );
// 	renderer.domElement.addEventListener( 'mousewheel', onDocumentMouseWheel, false );

// 	onMouseDownPosition = new THREE.Vector2();

// 	renderer.setClearColor( 0xffffff, 1);

// 	renderer.render(scene,camera);

// 	function updateCamera() {
// 		camera.position.x = camRadius * Math.sin( camTheta ) * Math.cos( camPhi );
// 		camera.position.y = camRadius * Math.sin( camPhi );
// 		camera.position.z = camRadius * Math.cos( camTheta ) * Math.cos( camPhi );
// 		camera.lookAt(new THREE.Vector3(0,1,0));
// 	}

// 	/* got the camera rotation code from: http://www.mrdoob.com/projects/voxels/#A/ */
// 	function onDocumentMouseDown( event ) {
// 		isMouseDown = true;
// 		onMouseDownTheta = camTheta;
// 		onMouseDownPhi = camPhi;
// 		onMouseDownPosition.x = event.clientX;
// 		onMouseDownPosition.y = event.clientY;
// 	}

// 	function onDocumentMouseMove( event ) {
// 		event.preventDefault();
// 		if ( isMouseDown ) {
// 			camTheta = - ( ( event.clientX - onMouseDownPosition.x ) * 0.01 ) + onMouseDownTheta;
			
// 			var dy = event.clientY - onMouseDownPosition.y;
			
// 			var newCamPhi = ( ( dy ) * 0.01 ) + onMouseDownPhi;

// 			if (newCamPhi < Math.PI/2 && newCamPhi > -Math.PI/2) {
// 				camPhi = newCamPhi;
// 			}
// 		}

// 		updateCamera();
// 		renderer.render(scene, camera);
// 	}

// 	function onDocumentMouseUp( event ) {
// 		event.preventDefault();
// 		isMouseDown = false;
// 	}

// 	function onDocumentMouseWheel( event ) { camRadius -= event.wheelDeltaY*.01; }

// 	function buildAxes( length ) {
// 	        var axes = new THREE.Object3D();

// 	        axes.util.add( buildAxis( new THREE.Vector3( 0, 0, 0 ), new THREE.Vector3( length, 0, 0 ), 0xFF0000 ) ); // +X
// 	        axes.util.add( buildAxis( new THREE.Vector3( 0, 0, 0 ), new THREE.Vector3( -length, 0, 0 ), 0xFF0000) ); // -X
// 	        axes.util.add( buildAxis( new THREE.Vector3( 0, 0, 0 ), new THREE.Vector3( 0, length, 0 ), 0x00FF00 ) ); // +Y
// 	        axes.util.add( buildAxis( new THREE.Vector3( 0, 0, 0 ), new THREE.Vector3( 0, -length, 0 ), 0x00FF00 ) ); // -Y
// 	        axes.util.add( buildAxis( new THREE.Vector3( 0, 0, 0 ), new THREE.Vector3( 0, 0, length ), 0x0000FF ) ); // +Z
// 	        axes.util.add( buildAxis( new THREE.Vector3( 0, 0, 0 ), new THREE.Vector3( 0, 0, -length ), 0x0000FF ) ); // -Z

// 	        return axes;

// 	}

// 	function buildAxis( src, dst, colorHex, dashed ) {
// 	        var geom = new THREE.Geometry(),
// 	            mat; 

// 	        if(dashed) {
// 	                mat = new THREE.LineDashedMaterial({ linewidth: 3, color: colorHex, dashSize: 3, gapSize: 3 });
// 	        } else {
// 	                mat = new THREE.LineBasicMaterial({ linewidth: 3, color: colorHex });
// 	        }

// 	        geom.vertices.push( src.clone() );
// 	        geom.vertices.push( dst.clone() );
// 	        geom.computeLineDistances(); // This one is SUPER important, otherwise dashed lines will appear as simple plain lines

// 	        var axis = new THREE.Line( geom, mat, THREE.LinePieces );

// 	        return axis;

// 	}

// 	this.animate = function(bouncePath,startTime) {
// 		if (startTime === undefined) {
// 			startTime = (new Date()).getTime();
// 		}

// 		var animationSpeed = 1;

// 		/* find time in the pattern and translate that to a discrete step in the prop position arrays */
// 		var timeElapsed = ((new Date()).getTime() - startTime)*animationSpeed/1000;
// 		var t = timeElapsed % bouncePath.T; // need to *1000 b/c timeElapsed is in ms
// 		var step = Math.floor(t/bouncePath.T*bouncePath.path.length);

// 		ballMesh.position = bouncePath.path[step];

// 		renderer.render(scene, camera);

// 		if (timeElapsed < bouncePath.T) {
// 			var self = this;
// 			requestAnimationFrame(function() { self.animate(bouncePath,startTime); });
// 		}
// 	}

// }

// function reportStats(ga) {
// 	var bestFitness;
// 	var T;
// 	if (ga.fittestMembers[ga.generations-1].fitness < 100) {
// 		bestFitness = ga.fittestMembers[ga.generations-1].fitness;
// 		T = ga.fittestMembers[ga.generations-1].T;
// 				var msg = 'Generation ' + ga.generations + ', best fitness: ' + bestFitness.toFixed(2) + ', ' + T.toFixed(2) + 's, velocity: ' + ga.fittestMembers[ga.generations-1].v.x.toFixed(2) + ', ' + ga.fittestMembers[ga.generations-1].v.y.toFixed(2) + ', ' + ga.fittestMembers[ga.generations-1].v.z.toFixed(2);
// 		$('#stats').append('<a href="#" onclick="animator.animate(ga.getBouncePath(ga.fittestMembers[' + (ga.generations-1).toString() + '].v))">'+msg+'</span><br/>');
// 	}

// 	if (ga.ableToFindSolution === undefined) {
// 		setTimeout(function() { reportStats(ga); },100);
// 	}	
// }

// function getConfigFromInput() {
// 	var input = $('#fitness').val().split('\n');
// 	/*
// 	p0.x,p0.y,p0.z,pT.x,pT.y,pT.z,T
// 	R,C
// 	numBounces,bounceOrder
// 	tossSign,catchSign
// 	surface1
// 	surface2
// 	surfaceN
// 	*/

// 	var posAndTime = input[0].split(',');
// 	var p0 = {x:parseFloat(posAndTime[0]), y:parseFloat(posAndTime[1]), z:parseFloat(posAndTime[2])};
// 	var pT = {x:parseFloat(posAndTime[3]), y:parseFloat(posAndTime[4]), z:parseFloat(posAndTime[5])};
// 	var minT = parseFloat(posAndTime[6]);
// 	var maxT = parseFloat(posAndTime[7]);

// 	var ballProperties = input[1].split(',');
// 	var R = parseFloat(ballProperties[0]);
// 	var C = parseFloat(ballProperties[1]);

// 	var bounces = input[2].split(',');
// 	var numBounces = parseInt(bounces[0]);
// 	var surfaceBounceOrder = undefined;
// 	if (bounces.length > 1) {
// 		surfaceBounceOrder = bounces.splice(1).map(function(a) { return parseInt(a); });
// 	}

// 	var tossCatchSigns = input[3].split(',');
// 	var tossSign = tossCatchSigns[0].trim() == 'u' ? undefined : parseInt(tossCatchSigns[0]);
// 	var catchSign = tossCatchSigns[1].trim() == 'u' ? undefined : parseInt(tossCatchSigns[1]);

// 	var surfaces = [];
// 	for (var i = 4; i < input.length; i++) {
// 		var surfaceProperties = input[i].split(',');
// 		var normal = {x:parseFloat(surfaceProperties[0]), y:parseFloat(surfaceProperties[1]), z:parseFloat(surfaceProperties[2])};
// 		var position = {x:parseFloat(surfaceProperties[3]), y:parseFloat(surfaceProperties[4]), z:parseFloat(surfaceProperties[5])};
// 		var scale = parseFloat(surfaceProperties[6]);
// 		surfaces.push({
// 			normal: normal,
// 			position: position,
// 			scale: scale
// 		});
// 	}

// 	var fitnessConfig = {
// 		p0: p0,
// 		pT: pT,
// 		minT: minT,
// 		maxT: maxT,
// 		R: R,
// 		C: C,
// 		dt: .01,
// 		G: -9.8,
// 		numBounces: numBounces,
// 		surfaceBounceOrder: surfaceBounceOrder,
// 		tossSign: tossSign,
// 		catchSign: catchSign,
// 		surfaces: surfaces
// 	};
	
// 	input = $('#ga').val().split(',');

// 	var gaConfig = {
// 		maxGenerations: parseInt(input[0]),
// 		populationSize: parseInt(input[1]),
// 		mutationChance: parseFloat(input[2]),
// 		mutationScale: parseFloat(input[3]),
// 		initialScale: parseFloat(input[4]),
// 		fitnessThreshold: parseFloat(input[5]),
// 		noGA: input.length > 6 && input[6] == 1 ? true : false
// 	}

// 	return {fitnessConfig: fitnessConfig, gaConfig: gaConfig};
// }

// var ga;
// var animator;

// function go() {
// 	var configs = getConfigFromInput();
// 	ga = new BounceGA(configs.gaConfig,configs.fitnessConfig);
// 	var done = ga.evolve();

// 	animator = new Animator(ga);

// 	function runAnimation() {
// 		if (ga.ableToFindSolution == true) {
// 			animator.animate(ga.getBouncePath(ga.fittestMembers[ga.fittestMembers.length-1].v));
// 		} else {
// 			setTimeout(runAnimation,0);
// 		}
// 	}

// 	$('#stats').empty();
// 	setTimeout(function() { reportStats(ga); },0);
	
// 	setTimeout(runAnimation,0);
// }

// function setExample(i) {
// 	switch (i) {
// 		case 0:
// 			$('#fitness').val("-.5,1.5,0, .5,1.5,0, 1.5,2\n\
// .05, .97\n\
// 1\n\
// 1, -1\n\
// 0,1,0, 0,0,0, 1");
// 			break;
// 		case 1:
// 			$('#fitness').val("-.5,1.5,0, .5,1.5,0, 2,2.1\n\
// .05, .97\n\
// 2, 0,1\n\
// u, u\n\
// .4,1,-.3, -1,0,1, 1\n\
// -.4,1,-.3, 1,0,1, 1");
// 			break;
// 		case 2:
// 			$('#fitness').val("-.5,1.5,0, .5,1.5,0, .3,1\n\
// .05, .97\n\
// 2\n\
// u, u\n\
// 1,0,0, 1,2,0, 2\n\
// 1,0,0, -1,2,0, 2");
// 			break; 
// 		case 3:
// 			$('#fitness').val("-.5,1.5,0, .5,1.5,0, 2,3\n\
// .05, .97\n\
// 3, 1,2,0\n\
// u, u\n\
// 1,0,0, 1,2,0, 1\n\
// 1,0,0, -1,2,0, 1\n\
// 0,1,0, 0,0,0, 1");
// 			break; 
// 	}
// }

// setExample(0);
},{"./util.js":7}],3:[function(require,module,exports){
var BounceGA = require('./BounceGA.js');
var Bezier = require('./Bezier.js');
var util = require('./util.js');

/* calculates the sum of all throws in the siteswap. used to determine the number of props */
function sumThrows(str) {

	var total = 0;
	for (var i = 0; i < str.length; i++) {
		if(parseInt(str[i])) {
			total += parseInt(str[i]);					
		} else if (str.charCodeAt(i) >= 97 && str.charCodeAt(i) <= 119) {
			// handle "a" through "z" (where "a" = 10)
			total += str.charCodeAt(i)-87;
		}

		// if the current character is a pass/spin marker
		// ignore the next character so we don't count the
		// juggler identifier  in something like <5p2|5p3|5p1>
		if ((str[i] == "P" || str[i] == "S") && parseInt(str[i+1]) ){
			i++;
		}
		// if the current character is a bounce marker
		// and then next character is a {, move forward until we find a }
		if ((str[i] == "B" || str[i] == "D" || str[i] == "T" || str[i] == "C" || str[i] == "S") && str[i+1] == "{") {
			i = str.indexOf("}",i);
		}
	}

	return total;
}

var flightPathCache = {};

/* CONSTANTS */

var LEFT = 0, RIGHT = 1;

/* core functions */
 function CreateSiteswap(siteswapStr, options) {
	
	/* return variable */
	var siteswap = {
		siteswap: 				siteswapStr,
		validSyntax: 			false,
		validPattern: 			false,
		collision:              undefined,
		multiplex: 				undefined,
		sync: 					undefined,
		pass: 					undefined,
		numJugglers: 			undefined,
		numProps: 				undefined,
		maxHeight: 				undefined,
		tosses: 				undefined,
		beats: 					undefined,
		states: 				undefined,
		propOrbits: 			undefined,
		propPositions: 			undefined,
		propRotations:  		undefined,
		jugglerHandPositions: 	undefined,
		jugglerElbowPositions: 	undefined,
		jugglers: 				undefined,
		validationOnly:			undefined,
		numStepsPerBeat:		undefined,		
		numSteps: 				undefined,
		beatDuration: 			undefined,
		dwellDuration: 			undefined,
		props:  				undefined,
		dwellPath: 				undefined,
		tossMatchVelocity:		undefined,
		catchMatchVelocity:		undefined,
		dwellCatchScale:		undefined,
		dwellTossScale:			undefined,
		emptyTossScale:			undefined,
		emptyCatchScale:		undefined,
		armAngle: 				undefined,
		surfaces: 				undefined,		
		errorMessage:  			undefined
	};

	/* regexps */
	var validTossRe,
		validMultiplexRe,
		validSyncRe,
		validBeatRe,
		validPassRe,
		validSiteswapRe;

	validateSyntax();

	setDefaultOptions();	
	
	if (siteswap.errorMessage) { return siteswap; }
	
	validatePattern();
	
	if (siteswap.errorMessage) { return siteswap; }

	if (!siteswap.validationOnly) {
		generatePropPositions();
	}

	return siteswap;

	function setDefaultOptions() {

		if (options === undefined) {
			options = {};
		}

		siteswap.validationOnly = (options.validationOnly === undefined ? false : options.validationOnly);		
		siteswap.beatDuration = (options.beatDuration === undefined ? .2 : options.beatDuration);		
		siteswap.dwellDuration = (options.dwellRatio === undefined ? siteswap.beatDuration*.5 : siteswap.beatDuration*options.dwellRatio);
		siteswap.numStepsPerBeat = (options.numStepsPerBeat === undefined ? Math.floor(siteswap.beatDuration*200) : options.numStepsPerBeat);
		siteswap.matchVelocity = (options.matchVelocity === undefined ? false : options.matchVelocity);
		siteswap.dwellCatchScale = (options.dwellCatchScale === undefined ? 0.05 : options.dwellCatchScale);
		siteswap.dwellTossScale = (options.dwellTossScale === undefined ? 0.05 : options.dwellTossScale);
		siteswap.emptyTossScale = (options.emptyTossScale === undefined ? 0.05 : options.emptyTossScale);
		siteswap.emptyCatchScale = (options.emptyCatchScale === undefined ? 0.05 : options.emptyCatchScale);
		siteswap.armAngle = (options.armAngle === undefined ? 0.1 : options.armAngle);
		
		if (options.props === undefined) {
			siteswap.props = [{type: 'ball', radius: .05, C: .95}];
		} else {
			siteswap.props = options.props;
		}

		if (options.dwellPath === undefined) {
			siteswap.dwellPath = [[{x:0.3,y:0,z:0,rotation:{x:4,y:0,z:-1,th:Math.PI/2}},{x:0.1,y:0,z:0,rotation:{x:4,y:0,z:1,th:Math.PI/2}}]];
		} else {
			siteswap.dwellPath = options.dwellPath;
		}

		if (options.surfaces === undefined) {
			siteswap.surfaces = [{position:{x:0,y:0,z:0}, normal:{x:0,y:1,z:0}, scale: 5}];
		} else {
			siteswap.surfaces = options.surfaces;
		}

		siteswap.jugglers = [];
		// if no juggler's were specified or there was a mismatch, just default to jugglers facing eachother
		if (options.jugglers === undefined || siteswap.numJugglers != options.jugglers.length) {				
			for (var i = 0; i < siteswap.numJugglers; i++) {
				siteswap.jugglers.push(
					{
						position: {x:0,z:-2*i},
						rotation: i*Math.PI
					}
				);
			}
		} else {
			for (var i = 0; i < options.jugglers.length; i++) {
				siteswap.jugglers.push(
					{
						position: options.jugglers[i].position,
						rotation: options.jugglers[i].rotation
					}
				);
			}
		}

		// set up axes on surfaces
		for (var i = 0; i < siteswap.surfaces.length; i++) {
			var surface = siteswap.surfaces[i];
			util.normalize(surface.normal);
			var axis1;
			if (surface.normal.x == 0 && surface.normal.z == 0) {
				axis1 = {x:1, y:0, z:0};
			} else {
				axis1 = {x:-surface.normal.z, y:0, z:surface.normal.x};
			}
			var axis2 = util.cross(surface.normal,axis1);
			util.normalize(axis1);
			util.multiply(axis1,surface.scale);
			util.normalize(axis2);
			util.multiply(axis2,surface.scale);
			surface.axis1 = axis1;
			surface.axis2 = axis2;
		}

	}

	/* check that the siteswap has the correct syntax */
	function validateSyntax() {
		var numJugglers = 1;
		var isPassingPattern = /<[^ ]+>/.test(siteswapStr);

		if (isPassingPattern) {
			var passingBeatArray = siteswapStr.match(/<[^ <>]+>/g);
			numJugglers = passingBeatArray[0].split("|").length;

			/* 
				check to make sure each beat in the passing pattern has the same number of jugglers 
				if a passing pattern only has 1 juggler than it's automatically a mismatch
			*/
			if(numJugglers == 1) {
				return siteswap;
			};
			
			var numJugglersTmp = numJugglers;
			passingBeatArray.map(function(a) { 
				if (a.split("|").length != numJugglersTmp) 
					{ return siteswap; } 
			});
		}

		/* the number of jugglers determines a valid pass pattern */
		var passPattern = "";
		if (numJugglers == 2) {
			passPattern = "P";
		} else if (numJugglers > 2) {
			passPattern = "P[1-" + numJugglers + "]";
		}

		/* construct the various regex patterns. see blog post for details about this */
		var validToss = "(R|L)?([\\da-o])x?A?(" + passPattern + ")?(C{(C|P)?})?(T{(C|P)?})?(B({\\d*(L|HL|F|HF)?})?)?(S{-?\\d+(.\\d+)?(,-?\\d+(.\\d+)?,-?\\d+(.\\d+)?,-?\\d+(.\\d+)?)?})?(D{\\d*\\.?\\d*})?";
		var validMultiplex = "\\[(" + validToss + ")+\\]";
		var validSync = "\\((" + validToss + "|" + validMultiplex + "),(" + validToss + "|" + validMultiplex + ")\\)";
		var validBeat = "(" + validToss + "|" + validMultiplex + "|" + validSync + ")";
		var validPass = "<" + validBeat + "(\\|" + validBeat + ")+>";
		var validSiteswap = "^(" + validPass + ")+|(" + validBeat + ")+$";

		validTossRe = new RegExp(validToss,"g");
		validMultiplexRe = new RegExp(validMultiplex,"g");
		validSyncRe = new RegExp(validSync,"g");
		validBeatRe = new RegExp(validBeat,"g");
		validPassRe = new RegExp(validPass,"g");
		validSiteswapRe = new RegExp(validSiteswap,"g");

		if (siteswapStr.match(validSiteswapRe) == siteswapStr) {
			siteswap.validSyntax = true;
			siteswap.multiplex = siteswapStr.match(validMultiplexRe) ? true : false;
			siteswap.sync = siteswapStr.match(validSyncRe) ? true : false;
			siteswap.pass = siteswapStr.match(validPassRe) ? true : false;
			siteswap.numJugglers = numJugglers;
		} else {
			siteswap.errorMessage = "Invalid syntax";
		} 
	}

	/* helper to get all the tosses for a given beat's siteswap */
	function getTosses(tosses, siteswapStr, juggler, sync, hand, dwellPathIx) {
		if (siteswapStr.match(validPassRe)) {
			var patterns = siteswapStr.match(validBeatRe);
			patterns.map(function(s,ix) {
				dwellPathIx = getTosses(tosses, s, ix, undefined, undefined, dwellPathIx);
				});
		} else if (siteswapStr.match(validSyncRe)) {
			var patterns = siteswapStr.split(",");
			dwellPathIx = getTosses(tosses,patterns[0].substr(1),juggler,true,LEFT, dwellPathIx);
			dwellPathIx = getTosses(tosses,patterns[1].substr(0,patterns[1].length-1),juggler,true,RIGHT, dwellPathIx);
		} else if (siteswapStr.match(validMultiplexRe)) {
			var patterns = siteswapStr.match(validTossRe);
			patterns.map(function(s,ix,arr) {
				dwellPathIx = getTosses(tosses,s,juggler, undefined, hand, dwellPathIx);
				if (ix < arr.length-1) {
					if (dwellPathIx == 0) {
						dwellPathIx = siteswap.dwellPath.length-1;
					} else {
						dwellPathIx--;
					}					
				}
			});
		} else {
			/* will work from "a" to "z" */
			var numBeats = (siteswapStr[0].charCodeAt(0) >= 97 && siteswapStr[0].charCodeAt(0) <= 119) ? siteswapStr[0].charCodeAt(0)-87 : parseInt(siteswapStr[0]);
			var targetJuggler = juggler;

			var pIx = siteswapStr.indexOf("P");
			var isPass = false;
			if (
				pIx > 0 &&
				siteswapStr[pIx+1] != "}" // check that the next character isn't a }, in which case this is a catch/toss penguin modifier
			) {				
				if (siteswap.numJugglers > 2) {					
					targetJuggler = parseInt(siteswapStr[pIx+1])-1;
				} else {
					targetJuggler = 1 - juggler;
				}
				isPass = true;
			}

			var numBounces = 0;
			var bounceOrder = [];
			var bIx = siteswapStr.indexOf("B");
			if (bIx > 0) {
				if (siteswapStr[bIx+1] == "{" && !isNaN(siteswapStr[bIx+2])) {
					numBounces = parseInt(siteswapStr[bIx+2]);
					for (var i = bIx + 3; i < siteswapStr.length; i++) {
						if (!isNaN(siteswapStr[i])) {
							var surfaceIx = parseInt(siteswapStr[i]);
							if (surfaceIx >= siteswap.surfaces.length) {
								throw {message: "Bounce surface index out of range"};
							} else {
								bounceOrder.push(surfaceIx);
							}							
						} else {
							break;
						}
					}
				} else {
					numBounces = 1;
					bounceOrder = [0];
				}
				if (bounceOrder.length < numBounces) {
					var numMissingBounces = bounceOrder.length;
					for (var i = 0; i < numBounces-numMissingBounces; i++) {
						bounceOrder.push(0);
					}
				}
			}

			var bounceType = undefined;
			if (numBounces > 0) {
				if (siteswapStr.match("HF")) {
					bounceType = "HF";
				} else if (siteswapStr.match("HL")) {
					bounceType = "HL";
				} else if (siteswapStr.match("F")) {
					bounceType = "F";
				} else if (siteswapStr.match("L")) {
					bounceType = "L";
				} else {
					bounceType = "L";
				}
			}

			var dIx = siteswapStr.indexOf("D");
			var dwellDuration;
			if (dIx > 0) {
				dwellDuration = siteswap.beatDuration*parseFloat(siteswapStr.substring(dIx+2,siteswapStr.indexOf("}")));
			}

			var tIx = siteswapStr.indexOf("T");
			var tossType = 'standard';
			if (tIx > 0) {
				var tossTypeId = siteswapStr.substring(tIx+2,siteswapStr.indexOf('}',tIx));
				if (tossTypeId.match("C")) {
					tossType = 'claw';
				} else if (tossTypeId.match("P")) {
					tossTypeId = "penguin";
				}
			}

			var cIx = siteswapStr.indexOf("C");
			var catchType = 'standard';
			if (cIx > 0) {
				var catchTypeId = siteswapStr.substring(cIx+2,siteswapStr.indexOf('}',cIx));
				if (catchTypeId.match("C")) {
					catchType = 'claw';
				} else if (catchTypeId.match("P")) {
					catchType = 'penguin';
				}
			}

			var crossing = numBeats % 2 == 1 ? true : false;
			// if the second character is an "x" then crossing is flipped
			if (siteswapStr.length > 1 && siteswapStr[1] == "x") {
				crossing = !crossing;
			}

			if (siteswapStr[0] == "R") {
				hand = RIGHT;
			} else if (siteswapStr[0] == "L") {
				hand = LEFT;
			}

			var numSpins;
			var tossOrientation = util.normalize({x:.1,y:.1,z:1});

			var sIx = siteswapStr.indexOf("S");			
			if (sIx > 0) {
				
				var spinConfig = siteswapStr.substring(sIx+2,siteswapStr.indexOf('}',sIx)).match(/-?\d+(\.\d+)?/g);				
				numSpins = parseFloat(spinConfig[0]);

				if (spinConfig.length > 1) {
					tossOrientation.x = parseFloat(spinConfig[1]);
					tossOrientation.y = parseFloat(spinConfig[2]);
					tossOrientation.z = parseFloat(spinConfig[3]);
					util.normalize(tossOrientation);
				}

			} else {
				numSpins = Math.floor(numBeats/2) + .2;
				// passes get an extra bit of spin
				if (isPass) {
					numSpins += .1;
				}

			}

			tosses.push(
				{
					juggler: juggler,
					targetJuggler: targetJuggler,
					hand: hand,
					crossing: crossing,
					numBeats: numBeats,
					siteswapStr: siteswapStr,
					numBounces: numBounces,
					bounceOrder: bounceOrder,
					bounceType: bounceType,
					numSpins: numSpins,					
					dwellPathIx: dwellPathIx,
					dwellDuration: dwellDuration === undefined ? siteswap.dwellDuration : dwellDuration,
					tossType: tossType,
					catchType: catchType,
					tossOrientation: tossOrientation,
					rotationAxis: {x:1,y:0,z:0},
					hold: numBeats == 2 && !crossing && siteswapStr.indexOf("A") == -1 ? true : false
				}
			);

			if (dwellPathIx == siteswap.dwellPath.length-1) {
				dwellPathIx = 0;
			} else {
				dwellPathIx++;
			}

		}		

		return dwellPathIx;
	}

	/* check that the siteswap is a repeatable pattern */
	function validatePattern() {

		/* get the array of each siteswap.beats' tosses */
		siteswap.beats = siteswap.pass ? siteswapStr.match(validPassRe) : siteswapStr.match(validBeatRe);		

		/* add (0,0) after each synchronous throw - this prevents the halving issue */
		for(var i = 0; i < siteswap.beats.length; i++) {
			if (siteswap.beats[i].match(validSyncRe)) {
				siteswap.beats.splice(i+1,0,'(0,0)');
				i++;
			}
		}
		
		/* figure out how many props */
		var tmp = 0;
		siteswap.beats.map(function(beat) {
			if (beat.match(validPassRe)) {
				var patterns = beat.split('|');
				for (var i = 0; i < patterns.length; i++) {
					if (i == 0) {
						patterns[i] = patterns[i].substr(1);
					} 
					if (i == patterns.length-1) {
						patterns[i] = patterns[i].substr(0,patterns[i].length-1);
					}
					tmp += sumThrows(patterns[i]);
				}
			} else {
				tmp += sumThrows(beat);
			}
		});

		if((tmp/siteswap.beats.length % 1) == 0 && tmp/siteswap.beats.length > 0) {
			siteswap.numProps = tmp/siteswap.beats.length;
		} else {		
			siteswap.errorMessage = "Cannot determine number of props";
			return;
		}

		/* make sure props array is correct length */
		while (siteswap.props.length < siteswap.numProps) {
			siteswap.props.push(util.cloneObject(siteswap.props.last()));
		}
		while (siteswap.props.length > siteswap.numProps) {
			siteswap.props.pop();
		}

		siteswap.tosses = [];

		var dwellPathIx = 0;

		/* for each beat get all the tosses */
		for (var i = 0; i < siteswap.beats.length; i++) {
			var tosses = [];
			dwellPathIx = getTosses(tosses,siteswap.beats[i], 0 /* assume juggler 0 */, undefined, undefined, dwellPathIx);
			siteswap.tosses.push(tosses);

			/* if the dwell paths aren't starting over at the same time as the beats, restart the pattern */
			if (i == siteswap.beats.length-1 && dwellPathIx != 0) {
				i = -1; // will get set to 0 above
			}

		}

		/* figure out the max throw height which will inform the size of the state array */
		siteswap.maxHeight = 0;

		for (var i = 0; i < siteswap.tosses.length; i++) {
			for (var j = 0; j < siteswap.tosses[i].length; j++) {
				if(siteswap.tosses[i][j].numBeats > siteswap.maxHeight) {
					siteswap.maxHeight = siteswap.tosses[i][j].numBeats;
				}
			}
		}

		/* ------------------------------------ */
		/* GENERATE STATE ARRAY AND PROP ORBITS */
		/* ------------------------------------ */

		/* create a queue of props */
		var props = [];

		for (var i = 0; i < siteswap.numProps; i++) {
			props.push(i);
		}

		/* initialize the state and prop orbits array */
		siteswap.states = [];
		siteswap.propOrbits = [];

		/* initialize current state */
		var curState = [];
		for (var j = 0; j < siteswap.numJugglers; j++) {
			curState.push([[],[]]);
			for (var k = 0; k < siteswap.maxHeight; k++) {
				curState[j][LEFT].push(undefined);
				curState[j][RIGHT].push(undefined);
			}
		}

		var patternComplete = false;
		var initComplete = false;
		var beat = 0;
		var hand = LEFT; /* default to starting with the left hand. this will automatically alternate */

		/* keep going until pattern complete */
		while (!patternComplete) {

			/* TODO: explain this */
			var tmpPropOrbits = util.cloneObject(siteswap.propOrbits);

			/* queue of props to throw this beat */
			var propsLanding = [];

			/* update the current state for each juggler */
			for (var j = 0; j < siteswap.numJugglers; j++) {
				var landingLeft = curState[j][LEFT].shift();
				if (landingLeft) {
					for (var k = 0; k < landingLeft.length; k++) {
 						propsLanding.push({propId: landingLeft[k], juggler: j, hand: LEFT});	
					}						
				}
				var landingRight = curState[j][RIGHT].shift();
				if (landingRight) {
					for (var k = 0; k < landingRight.length; k++) {
						propsLanding.push({propId: landingRight[k], juggler: j, hand: RIGHT});	
					}						
				}					
				curState[j][LEFT].push(undefined);
				curState[j][RIGHT].push(undefined);
			}

			/* iterate through all the tosses and update the current state */
			for (var j = 0; j < siteswap.tosses[beat % siteswap.tosses.length].length; j++) {
				
				var toss = siteswap.tosses[beat % siteswap.tosses.length][j];
				var tossHand = (toss.hand == undefined ? hand : toss.hand);
				var catchHand = (toss.crossing ? 1 - tossHand : tossHand);

				var prop = undefined;

				/* iterate through the props landing and look for one landing in the hand that this toss is occurring */
				for (var k = 0; k < propsLanding.length; k++) {
					if(propsLanding[k].juggler == toss.juggler && propsLanding[k].hand == tossHand) {
						
						/* if a prop is landing in a hand this is tossing a 0 then invalid siteswap */
						if (toss.numBeats == 0) {
							siteswap.errorMessage = "Prop landing on 0 toss at beat " + beat;
							return;
						}

						prop = propsLanding.splice(k,1)[0].propId;
						break;
					}
				}

				/* if no props landing to be thrown, get one from the queue - only if this isn't a 0 toss */
				if (prop == undefined && toss.numBeats > 0) {
					prop = props.shift();			
				} 

				/* if prop is still undefined (ie. there are none left) then we've got an invalid siteswap - only if this isn't a 0 toss */
				if (prop == undefined && toss.numBeats > 0) {
					siteswap.errorMessage = "No prop available to toss at beat " + beat;
					return;
				}

				/* so long as this isn't a 0 toss, update the current state and append to prop orbits */
				if (toss.numBeats > 0) {
					
					if(!tmpPropOrbits[prop]) {
						tmpPropOrbits[prop] = [];
					}

					tmpPropOrbits[prop].push({beat: beat, numBeats: toss.numBeats, juggler: toss.juggler, hand: tossHand, numBounces: toss.numBounces, bounceType: toss.bounceType, bounceOrder: toss.bounceOrder, numSpins: toss.numSpins, dwellPathIx: toss.dwellPathIx, dwellDuration: toss.dwellDuration, tossType: toss.tossType, catchType: toss.catchType, tossOrientation: toss.tossOrientation, rotationAxis: toss.rotationAxis, hold: toss.hold });

					if(curState[toss.targetJuggler][catchHand][toss.numBeats-1] == undefined) {
						curState[toss.targetJuggler][catchHand][toss.numBeats-1] = [prop];
					} else {
						curState[toss.targetJuggler][catchHand][toss.numBeats-1].push(prop);
					}

				}
				
			}
							

			/* if we're at the beginning of the toss array and we've returned to the original state, the pattern is complete */
			if (initComplete && beat % siteswap.tosses.length == 0 && util.arraysEqual(siteswap.states[0],curState)) {
				patternComplete = true;				
			} else {
				/* add the current state to the state array and update prop orbits */
				siteswap.states.push(util.cloneObject(curState));
				siteswap.propOrbits = tmpPropOrbits;
			}					

			/* if all props have been introduced to pattern and we're at the end of the pattern, init is complete and steady-state pattern truly begins with the next beat */
			if (props.length == 0 && (beat+1) % siteswap.tosses.length == 0 && !initComplete) {
				initComplete = true;
				beat = -1;
				siteswap.states = []; /* reset the states and prop orbits */
				siteswap.propOrbits = [];
			}			

			beat++;
			hand = 1 - hand; //alternate hands

			/* fail safe in case the pattern is too long */
			if (beat > 1000) {
				siteswap.errorMessage = "Pattern took more than 1000 beats to repeat states"
				return;
			}

		}

		/* if we've gotten to this point, the pattern is repeatable and thus valid */
		siteswap.numSteps = siteswap.states.length*siteswap.numStepsPerBeat;
		siteswap.validPattern = true;
	}

	function generatePropPositions() {

		try {

			// clear flight path cache
			flightPathCache = {};

			/* initialize prop positions */
			var propPositions = [];
			for (var i = 0; i < siteswap.numProps; i++) {
				propPositions.push([]);
			}

			/* init prop rotations */
			var propRotations = [];
			for (var i = 0; i < siteswap.numProps; i++) {
				propRotations.push([]);
			}

			/* initialize juggler hand positions */
			var jugglerHandPositions = [];
			for (var i = 0; i < siteswap.numJugglers; i++) {
				jugglerHandPositions.push([[],[]]);
			}
			var jugglerElbowPositions = [];
			for (var i = 0; i < siteswap.numJugglers; i++) {
				jugglerElbowPositions.push([[],[]]);
			}

			/* generate prop positions */
			for (var step = 0; step < siteswap.numSteps; step++) {
				
				var tmpJugglerHandPositions = [];
				for (var i = 0; i < siteswap.numJugglers; i++) {
					tmpJugglerHandPositions.push([undefined,undefined]);
				}

				var currentBeat = Math.floor(step*siteswap.states.length/siteswap.numSteps);
				var currentTime = siteswap.beatDuration*step*siteswap.states.length/siteswap.numSteps;

				/* find the current state of each prop */
				for(var prop = 0; prop < siteswap.numProps; prop++) {					

					var prevToss = undefined, curToss = undefined, nextToss = undefined;
					
					if (siteswap.propOrbits[prop].length == 1) {
						
						prevToss = siteswap.propOrbits[prop][0];
						curToss = siteswap.propOrbits[prop][0];
						nextToss = siteswap.propOrbits[prop][0];

					}
					var orbitBeatFound = false;
					for (var i = 0; i < siteswap.propOrbits[prop].length-1; i++) {
						if (!orbitBeatFound && siteswap.propOrbits[prop][i].beat <= currentBeat && siteswap.propOrbits[prop][i+1].beat > currentBeat) {
							
							if (i == 0) {
								prevToss = siteswap.propOrbits[prop][siteswap.propOrbits[prop].length-1];
							} else {
								prevToss = siteswap.propOrbits[prop][i-1];
							}
							curToss = siteswap.propOrbits[prop][i];
							nextToss = siteswap.propOrbits[prop][i+1];

							orbitBeatFound = true;

						} else if (!orbitBeatFound && i == siteswap.propOrbits[prop].length-2) { 

							prevToss = siteswap.propOrbits[prop][i];
							curToss = siteswap.propOrbits[prop][i+1];
							nextToss = siteswap.propOrbits[prop][0];

						}
					}

					var tossTime = curToss.beat*siteswap.beatDuration+curToss.dwellDuration;
					var catchTime = nextToss.beat*siteswap.beatDuration;
					if (tossTime > catchTime && catchTime <= currentTime) {
						catchTime += (siteswap.beatDuration*siteswap.states.length);	
					}					
					else if (tossTime > catchTime && catchTime > currentTime) { 
						tossTime -= (siteswap.beatDuration*siteswap.states.length);
					}

					var lastTossTime = prevToss.beat*siteswap.beatDuration+prevToss.dwellDuration;
					var lastCatchTime = curToss.beat*siteswap.beatDuration;
					if (lastTossTime > lastCatchTime && lastCatchTime <= currentTime) {
						lastCatchTime += (siteswap.beatDuration*siteswap.states.length);	
					}
					else if (lastTossTime > lastCatchTime && lastCatchTime > currentTime) { 
						lastTossTime -= (siteswap.beatDuration*siteswap.states.length);
					}

					if (currentTime < tossTime) {
						/* interpolate dwell path */

						var launches = [];
						var landings = [];

						// iterate over all other props to see if any others are in the hand now as well
						if (siteswap.multiplex) {
							for (var i = 0; i < siteswap.propOrbits.length; i++) {
								if (i != prop) {
									for (j = 0; j < siteswap.propOrbits[i].length; j++) {
										var multiplexCurToss = siteswap.propOrbits[i][j];

										if (curToss.beat == multiplexCurToss.beat && curToss.juggler == multiplexCurToss.juggler && curToss.hand == multiplexCurToss.hand) {
											
											var multiplexNextToss = j == siteswap.propOrbits[i].length-1 ? siteswap.propOrbits[i][0] : siteswap.propOrbits[i][j+1];
											var multiplexPrevToss = j == 0 ? siteswap.propOrbits[i][siteswap.propOrbits[i].length-1] : siteswap.propOrbits[i][j-1];

											var multiplexTossTime = multiplexCurToss.beat*siteswap.beatDuration+multiplexCurToss.dwellDuration;
											var multiplexCatchTime = multiplexNextToss.beat*siteswap.beatDuration;
											if (multiplexTossTime > multiplexCatchTime && multiplexCatchTime <= currentTime) {
												multiplexCatchTime += (siteswap.beatDuration*siteswap.states.length);	
											}					
											else if (multiplexTossTime > multiplexCatchTime && multiplexCatchTime > currentTime) { 
												multiplexTossTime -= (siteswap.beatDuration*siteswap.states.length);
											}

											var multiplexLastTossTime = multiplexPrevToss.beat*siteswap.beatDuration+multiplexPrevToss.dwellDuration;
											var multiplexLastCatchTime = multiplexCurToss.beat*siteswap.beatDuration;
											if (multiplexLastTossTime > multiplexLastCatchTime && multiplexLastCatchTime <= currentTime) {
												multiplexLastCatchTime += (siteswap.beatDuration*siteswap.states.length);	
											}
											else if (multiplexLastTossTime > multiplexLastCatchTime && multiplexLastCatchTime > currentTime) { 
												multiplexLastTossTime -= (siteswap.beatDuration*siteswap.states.length);
											}

											launches.push(interpolateFlightPath(
												getDwellPosition(siteswap.dwellPath[multiplexCurToss.dwellPathIx],multiplexCurToss.juggler,multiplexCurToss.hand,1), /* p0 */
												getDwellPosition(siteswap.dwellPath[multiplexNextToss.dwellPathIx],multiplexNextToss.juggler,multiplexNextToss.hand,0), /* p1 */
												(multiplexCatchTime - multiplexTossTime),
												0,								
												{
													numBounces: multiplexCurToss.numBounces, 
													bounceType: multiplexCurToss.bounceType, 
													bounceOrder: multiplexCurToss.bounceOrder,
													R: siteswap.props[i].radius, 
													C: siteswap.props[i].C
												}
											));

											landings.push(interpolateFlightPath(
												getDwellPosition(siteswap.dwellPath[multiplexPrevToss.dwellPathIx],multiplexPrevToss.juggler,multiplexPrevToss.hand,1), /* p0 */
												getDwellPosition(siteswap.dwellPath[multiplexCurToss.dwellPathIx],multiplexCurToss.juggler,multiplexCurToss.hand,0), /* p1 */
												(multiplexLastCatchTime - multiplexLastTossTime),
												(multiplexLastCatchTime - multiplexLastTossTime),
												{
													numBounces: multiplexPrevToss.numBounces, 
													bounceType: multiplexPrevToss.bounceType,
													bounceOrder: multiplexPrevToss.bounceOrder, 
													R: siteswap.props[i].radius, 
													C: siteswap.props[i].C
												}
											));										
										}
									}
								}
							}
						}						

						launches.push(interpolateFlightPath(
								getDwellPosition(siteswap.dwellPath[curToss.dwellPathIx],curToss.juggler,curToss.hand,1), /* p0 */
								getDwellPosition(siteswap.dwellPath[nextToss.dwellPathIx],nextToss.juggler,nextToss.hand,0), /* p1 */
								(catchTime - tossTime),
								0,								
								{
									numBounces: curToss.numBounces, 
									bounceType: curToss.bounceType, 
									bounceOrder: curToss.bounceOrder,
									R: siteswap.props[prop].radius, 
									C: siteswap.props[prop].C
								}
							));

						landings.push(interpolateFlightPath(
								getDwellPosition(siteswap.dwellPath[prevToss.dwellPathIx],prevToss.juggler,prevToss.hand,1), /* p0 */
								getDwellPosition(siteswap.dwellPath[curToss.dwellPathIx],curToss.juggler,curToss.hand,0), /* p1 */
								(lastCatchTime - lastTossTime),
								(lastCatchTime - lastTossTime),
								{
									numBounces: prevToss.numBounces, 
									bounceType: prevToss.bounceType,
									bounceOrder: prevToss.bounceOrder, 
									R: siteswap.props[prop].radius, 
									C: siteswap.props[prop].C
								}
							));

						var land = {dx: 0, dy: 0, dz: 0, x: 0, y: 0, z:0};
						var launch = {dx: 0, dy: 0, dz: 0, x: 0, y: 0, z:0};
						for (i = 0; i < landings.length; i++) {
							land.dx += landings[i].dx/landings.length;
							land.dy += landings[i].dy/landings.length;
							land.dz += landings[i].dz/landings.length;
							land.x += landings[i].x/landings.length;
							land.y += landings[i].y/landings.length;
							land.z += landings[i].z/landings.length;

							launch.dx += launches[i].dx/landings.length;
							launch.dy += launches[i].dy/landings.length;
							launch.dz += launches[i].dz/landings.length;
							launch.x += launches[i].x/landings.length;
							launch.y += launches[i].y/landings.length;
							launch.z += launches[i].z/landings.length;
						}

						//var land = landings.last();
						//var launch = launches.last();

						var t = 1-(tossTime - currentTime)/curToss.dwellDuration;
						var pos = getDwellPosition(
							siteswap.dwellPath[curToss.dwellPathIx]
							, curToss.juggler
							, curToss.hand
							, t
							, land
							, launch
							, siteswap.dwellCatchScale
							, siteswap.dwellTossScale
						);

						// the landing return by flight path interpolation may be slightly off (if solved by BounceGA)
						// in this case we should find the correct landing and interpolate between the two 
						var correctLand = getDwellPosition(siteswap.dwellPath[curToss.dwellPathIx],curToss.juggler,curToss.hand,0);

						var landingDiff = {x: land.x - correctLand.x, y: land.y - correctLand.y, z: land.z - correctLand.z};
						pos.x += (1-t)*landingDiff.x;
						pos.y += (1-t)*landingDiff.y;
						pos.z += (1-t)*landingDiff.z;
						var catchAngle = Math.atan2(-land.dx,-land.dy);
						var tossAngle = Math.atan2(launch.dx,launch.dy);
						if (curToss.tossType == 'claw') {
							tossAngle -= Math.PI;
						} else if (curToss.tossTypeId == 'penguin') {
							tossAngle -= 2*Math.PI;
						}
						if (curToss.catchType == 'claw') {
							catchAngle -= Math.PI;
						} else if (curToss.catchType == 'penguin') {
							catchAngle -= 2*Math.PI;
						}
						pos.angle = catchAngle + t*(tossAngle-catchAngle);						
						if (curToss.hand == RIGHT)
							pos.angle *= -1;

						pos.dwell = true;
						
						propPositions[prop].push(pos);
						
						/* assign juggler hand positions */
						if (tmpJugglerHandPositions[curToss.juggler][curToss.hand] == undefined) {
							tmpJugglerHandPositions[curToss.juggler][curToss.hand] = pos;
						}					

						var q = getPropQuaternion(prevToss.tossOrientation, prevToss.rotationAxis, siteswap.jugglers[prevToss.juggler].rotation, prevToss.numSpins*2*Math.PI, prevToss.hand);
						var q2 = getPropQuaternion(curToss.tossOrientation, curToss.rotationAxis, siteswap.jugglers[curToss.juggler].rotation, 0, curToss.hand);
						q.slerp(q2, t);
						propRotations[prop].push(q);
					} else {

						/*
						calculate position at current time
						*/

						var T = catchTime - tossTime;
						var t = currentTime - tossTime;						
						var pos;

						// if the current toss is held (ie. a 2) then don't leave the hand 
						// so this is the same code as the empty hand's path
						if (curToss.hold) {
							//getDwellPosition(P,juggler,hand,t,v_0,v_T,v_0scale,v_Tscale,matchVelocity)		

							var p0 = getDwellPosition(siteswap.dwellPath[curToss.dwellPathIx],curToss.juggler,curToss.hand,1);
							var pT = getDwellPosition(siteswap.dwellPath[nextToss.dwellPathIx],nextToss.juggler,nextToss.hand,0);

							var v_0 = interpolateFlightPath(
								p0,
								pT,
								T,
								0
							);
							var v_T = interpolateFlightPath(
								p0,
								pT,
								T,
								T
							);
							
							pos = getDwellPosition(
								[siteswap.dwellPath[curToss.dwellPathIx].last(),siteswap.dwellPath[nextToss.dwellPathIx][0]]
								, curToss.juggler
								, curToss.hand
								, t/T
								, v_0
								, v_T
								, siteswap.emptyTossScale
								, siteswap.emptyCatchScale
							);

						} else {

							// if not holding prop then interpolate flight path

							pos = interpolateFlightPath(
								getDwellPosition(siteswap.dwellPath[curToss.dwellPathIx],curToss.juggler,curToss.hand,1), /* p0 */
								getDwellPosition(siteswap.dwellPath[nextToss.dwellPathIx],nextToss.juggler,nextToss.hand,0), /* p1 */
								T,
								t,
								{
									numBounces: curToss.numBounces, 
									bounceType: curToss.bounceType, 
									bounceOrder: curToss.bounceOrder,
									R: siteswap.props[prop].radius, 
									C: siteswap.props[prop].C
								}
							);
						}						

						pos.dwell = false;

						propPositions[prop].push(pos);

						var catchRotation = curToss.numSpins*2*Math.PI;
						var tossRotation = 0;
						var currentRotation = tossRotation + (t/T)*(catchRotation - tossRotation);

						propRotations[prop].push(getPropQuaternion(curToss.tossOrientation, curToss.rotationAxis, siteswap.jugglers[curToss.juggler].rotation, currentRotation, curToss.hand));

					}					

				}

				/* set hand positions that weren't set */
				for (var juggler = 0; juggler < siteswap.numJugglers; juggler++) {
					for (var hand = 0; hand <= 1; hand++) {
						if(tmpJugglerHandPositions[juggler][hand] == undefined) {
							
							/* need 
								nextToss - to determine where the hand is going to
								propLastToss - to determine where the prop we're catching came from so we know its catch velocity
								lastToss - to determine where the hand is coming from
								propNextToss - to determine where the prop we just tossed is going so we know its toss velocity
							*/


							/* find the next beat a prop is going to be in this hand and linearly move to the catch position */
							var nextToss = undefined, minToss = undefined, lastToss = undefined, maxToss = undefined;
							var minTossProp = undefined, minTossOrbit = undefined, nextTossProp = undefined, nextTossOrbit = undefined, maxTossProp = undefined, maxTossOrbit = undefined, lastTossProp = undefined, lastTossOrbit = undefined;
							var propLastToss = undefined, propNextToss = undefined;

							for (var prop = 0; prop < siteswap.propOrbits.length; prop++) {
								for (var orbit = 0; orbit < siteswap.propOrbits[prop].length; orbit++) {
									if (siteswap.propOrbits[prop][orbit].juggler == juggler && siteswap.propOrbits[prop][orbit].hand == hand) {
										
										// min beat
										if (minToss == undefined || siteswap.propOrbits[prop][orbit].beat < minToss.beat) {
											minToss = siteswap.propOrbits[prop][orbit];
											minTossProp = prop;
											minTossOrbit = orbit;
										}
										// next beat
										if (siteswap.propOrbits[prop][orbit].beat > currentBeat && (nextToss == undefined || siteswap.propOrbits[prop][orbit].beat < nextToss.beat)) {
											nextToss = siteswap.propOrbits[prop][orbit];
											nextTossProp = prop;
											nextTossOrbit = orbit;
										}
										
										// max beat
										if (maxToss == undefined || siteswap.propOrbits[prop][orbit].beat > maxToss.beat) {
											maxToss = siteswap.propOrbits[prop][orbit];
											maxTossProp = prop;
											maxTossOrbit = orbit;
										}
										// last beat
										if (siteswap.propOrbits[prop][orbit].beat <= currentBeat && (lastToss == undefined || siteswap.propOrbits[prop][orbit].beat > lastToss.beat)) {
											lastToss = siteswap.propOrbits[prop][orbit];
											lastTossProp = prop;
											lastTossOrbit = orbit;
										}
									}
								}
							}

							if (nextToss == undefined) {
								nextToss = minToss;
								nextTossProp = minTossProp;
								nextTossOrbit = minTossOrbit;
							}
							if (lastToss == undefined) {
								lastToss = maxToss;
								lastTossProp = maxTossProp;
								lastTossOrbit = maxTossOrbit;
							}

							if (nextTossOrbit == 0) {
								propLastToss = siteswap.propOrbits[nextTossProp].last();
							} else {
								propLastToss = siteswap.propOrbits[nextTossProp][nextTossOrbit-1];
							}

							if (lastTossOrbit == siteswap.propOrbits[lastTossProp].length-1) {
								propNextToss = siteswap.propOrbits[lastTossProp][0];
							} else {
								propNextToss = siteswap.propOrbits[lastTossProp][lastTossOrbit+1];
							}

							var nextCatchTime = nextToss.beat*siteswap.beatDuration;
							if (nextCatchTime < currentTime) {
								nextCatchTime += (siteswap.beatDuration*siteswap.states.length);
							}
							var lastThrowTime = lastToss.beat*siteswap.beatDuration+lastToss.dwellDuration;
							if (lastThrowTime > currentTime) {
								lastThrowTime -= (siteswap.beatDuration*siteswap.states.length);
							}
							var propNextCatchTime = propNextToss.beat*siteswap.beatDuration;
							if (propNextCatchTime < lastThrowTime) {
								propNextCatchTime += (siteswap.beatDuration*siteswap.states.length);
							}
							var propLastThrowTime = propLastToss.beat*siteswap.beatDuration+propLastToss.dwellDuration;
							if (propLastThrowTime > nextCatchTime) {
								propLastThrowTime -= (siteswap.beatDuration*siteswap.states.length);
							}

							var v_0 = interpolateFlightPath(
								getDwellPosition(siteswap.dwellPath[lastToss.dwellPathIx],lastToss.juggler,lastToss.hand,1), /* p0 */
								getDwellPosition(siteswap.dwellPath[propNextToss.dwellPathIx],propNextToss.juggler,propNextToss.hand,0), /* p1 */
								lastToss.numBeats*siteswap.beatDuration-lastToss.dwellDuration,
								0,
								{
									numBounces: lastToss.numBounces, 
									bounceType: lastToss.bounceType,
									bounceOrder: lastToss.bounceOrder,
									R: siteswap.props[lastTossProp].radius, 
									C: siteswap.props[lastTossProp].C
								}
							);
							var v_T = interpolateFlightPath(
								getDwellPosition(siteswap.dwellPath[propLastToss.dwellPathIx],propLastToss.juggler,propLastToss.hand,1), /* p0 */
								getDwellPosition(siteswap.dwellPath[nextToss.dwellPathIx],nextToss.juggler,nextToss.hand,0), /* p1 */
								propLastToss.numBeats*siteswap.beatDuration-propLastToss.dwellDuration,
								propLastToss.numBeats*siteswap.beatDuration-propLastToss.dwellDuration,
								{
									numBounces: propLastToss.numBounces, 
									bounceType: propLastToss.bounceType,
									bounceOrder: propLastToss.bounceOrder,
									R: siteswap.props[nextTossProp].radius, 
									C: siteswap.props[nextTossProp].C
								}
							);

							var t = (currentTime - lastThrowTime)/(nextCatchTime - lastThrowTime);
							var pos = getDwellPosition(
								[siteswap.dwellPath[lastToss.dwellPathIx].last(),siteswap.dwellPath[nextToss.dwellPathIx][0]]
								, lastToss.juggler
								, lastToss.hand
								, t
								, v_0
								, v_T
								, siteswap.emptyTossScale
								, siteswap.emptyCatchScale
							);

							var correctCatch = getDwellPosition(
								[siteswap.dwellPath[lastToss.dwellPathIx].last(),siteswap.dwellPath[nextToss.dwellPathIx][0]]
								, lastToss.juggler
								, lastToss.hand
								, 1
							);							

							var catchDiff = {x: v_T.x - correctCatch.x, y: v_T.y - correctCatch.y, z: v_T.z - correctCatch.z};
							pos.x += (t)*catchDiff.x;
							pos.y += (t)*catchDiff.y;
							pos.z += (t)*catchDiff.z;						

							var catchAngle = Math.atan2(-v_T.dx,-v_T.dy);
							var tossAngle = Math.atan2(v_0.dx,v_0.dy);
							if (lastToss.tossType == 'claw') {
								tossAngle -= Math.PI;
							} else if (lastToss.tossType == 'penguin') {
								tossAngle -= 2*Math.PI;
							}
							if (nextToss.catchType == 'claw') {
								catchAngle -= Math.PI;
							} else if (nextToss.catchType == 'penguin') {
								catchAngle -= 2*Math.PI;
							}
							pos.angle = tossAngle + t*(catchAngle-tossAngle);
							if (hand == RIGHT) {
								pos.angle *= -1;
							}

							tmpJugglerHandPositions[juggler][hand] = pos;
						}					

						jugglerHandPositions[juggler][hand].push(tmpJugglerHandPositions[juggler][hand]);
						jugglerElbowPositions[juggler][hand].push(
							getElbowPosition(
								{x:siteswap.jugglers[juggler].position.x+Math.cos(siteswap.jugglers[juggler].rotation)*(hand == LEFT ? - 1 : 1)*.225,y:1.425,z:siteswap.jugglers[juggler].position.z+Math.sin(siteswap.jugglers[juggler].rotation)*(hand == LEFT ? - 1 : 1)*.225}, // shoulder
								tmpJugglerHandPositions[juggler][hand], // hand position
								.45, // half arm length
								siteswap.armAngle, // chicken wing factor
								hand // hand
							)
						);						

					}

				}

			}

			siteswap.propPositions = propPositions;
			siteswap.propRotations = propRotations;
			siteswap.jugglerHandPositions = jugglerHandPositions;
			siteswap.jugglerElbowPositions = jugglerElbowPositions;

			siteswap.collision = checkForCollision();

		} catch(e) {
			siteswap.errorMessage = e.message;
		}
	}

	/* interpolate flight path */	
	function interpolateFlightPath(p0, p1, T, t, options) {
		/*
		p0 - starting position
		p1 - ending position
		T - total time
		t - time elapsed
		options - configurable parameters
			- bounceType - L,F,HL,HF
			- dt - time increment for bounce simulation
			- eps - margin of error acceptable for bouncing
			- dv - velocity increment for bouncing
			- numBounces
			- G - gravity
			- C - coefficient of restitution
			- tries - # of times to try and get it right before throwing error
			- R - radius of prop
		*/

		// round T to 2 decimal places for the sake of the flight path cache
		T = parseFloat(T.toFixed(2));

		if (options == undefined) { options = {}; }

		/* set defaults */
		if (!options.bounceType)
			options.bounceType = "L";
		if (!options.dt)
			options.dt = .01;
		if (!options.eps) 
			options.eps = .01;
		if (!options.dv)
			options.dv = .01;
		if (!options.numBounces)
			options.numBounces = 0;
		if (!options.bounceOrder)
			options.bounceOrder = [0];
		if (!options.G)
			options.G = -9.8;
		if (!options.C)
			options.C = .95;
		if (!options.tries) 
			options.tries = 10000;
		if (!options.R) /* should this one be required? */
			options.R = .1;

		var inputKey = JSON.stringify({p0:p0,p1:p1,T:T,options:options});

		if (options.numBounces == 0) {
			
			return  {
				x: p0.x + (p1.x-p0.x)*t/T,
				y: p0.y + (p1.y - p0.y - .5*options.G*T*T)*t/T + .5*options.G*t*t,
				z: p0.z + (p1.z-p0.z)*t/T,
				dx: (p1.x-p0.x)/T,
				dy: (p1.y - p0.y -.5*options.G*T*T)/T + options.G*t,
				dz: (p1.z-p0.z)/T
			};

		} else if (flightPathCache[inputKey] == undefined) {

			var fitnessConfig = {
				p0: p0,
				pT: p1,
				minT: T,
				maxT: T,
				R: options.R,
				C: options.C,
				dt: options.dt,
				G: -9.8,
				numBounces: options.numBounces,
				surfaceBounceOrder: options.bounceOrder,
				tossSign: options.bounceType == 'L' || options.bounceType == 'HL' ? 1 : -1,
				catchSign: options.bounceType == 'L' || options.bounceType == 'HF' ? 1 : -1,
				surfaces: siteswap.surfaces
			};

			var gaConfig = {
				maxGenerations: 500,
				populationSize: 50,
				mutationChance: .7,
				mutationScale: 5,
				initialScale: 20,
				fitnessThreshold: .05,
				noGA: false
			};

			ga = new BounceGA.BounceGA(gaConfig,fitnessConfig);
			ga.evolve();

			if (!ga.ableToFindSolution) {
				/* TODO - improve error to explain why the bounce path couldn't be calculated */
				throw {message: 'Unable to calculate bounce path'};
			}

			var gaResult = ga.getBouncePath(ga.fittestMembers[ga.fittestMembers.length-1].v);

			flightPathCache[inputKey] = {path: gaResult.path, velocities: gaResult.velocities};

		}

		var flightPath = flightPathCache[inputKey];
		return {
			x: flightPath.path[Math.floor((flightPath.path.length-1)*t/T)].x,
			y: flightPath.path[Math.floor((flightPath.path.length-1)*t/T)].y,
			z: flightPath.path[Math.floor((flightPath.path.length-1)*t/T)].z,
			dx: flightPath.velocities[Math.floor((flightPath.velocities.length-1)*t/T)].x,
			dy: flightPath.velocities[Math.floor((flightPath.velocities.length-1)*t/T)].y,
			dz: flightPath.velocities[Math.floor((flightPath.velocities.length-1)*t/T)].z,
		};
	}

	function getDwellPosition(P,juggler,hand,t,v_0,v_T,v_0scale,v_Tscale,matchVelocity) {

		if (hand == LEFT && (v_0 || v_T)) {
			v_0.dx *= -1;
			v_T.dx *= -1;
		}

		var dwellPosition = Bezier.interpolateBezierSpline(P,t,v_0,v_T,v_0scale,v_Tscale,matchVelocity);

		return {
			x: siteswap.jugglers[juggler].position.x + ((hand == LEFT ? -1 : 1)*dwellPosition.x)*Math.cos(siteswap.jugglers[juggler].rotation) - (dwellPosition.z - .4125)*Math.sin(siteswap.jugglers[juggler].rotation),
			y: 1.0125 + dwellPosition.y,
			z: siteswap.jugglers[juggler].position.z + ((hand == LEFT ? -1 : 1)*dwellPosition.x)*Math.sin(siteswap.jugglers[juggler].rotation) + (dwellPosition.z - .4125)*Math.cos(siteswap.jugglers[juggler].rotation)
		};
	}

	function getElbowPosition(S,H,l,w,hand) {
		var Hp = {};
		Hp.x = H.x - S.x;
		Hp.y = H.y - S.y;
		Hp.z = H.z - S.z;

		var Hpp = {};
		Hpp.x = Math.sqrt(Hp.x*Hp.x + Hp.z*Hp.z);
		Hpp.y = Hp.y;
		Hpp.z = 0;

		var th = Math.atan2(Hp.z,Hp.x);

		var magHp = Math.sqrt(Hp.x*Hp.x + Hp.y*Hp.y + Hp.z*Hp.z);

		/* magically stretch arms */
		if (2*l < magHp) {
			l = magHp/2;
		}

		var u1 = {};
		u1.x = Hpp.y/magHp;
		u1.y = -Hpp.x/magHp;
		u1.z = 0;

		var u2 = {x:0,y:0};
		if (hand == 1) {
			u2.z = -1;
		} else {
			u2.z = 1;
		}

		var h = Math.sqrt(l*l - .25*magHp*magHp);

		var Epp = {};
		Epp.x = Hpp.x/2 + h*u1.x*Math.cos(w) + h*u2.x*Math.sin(w);
		Epp.y = Hpp.y/2 + h*u1.y*Math.cos(w) + h*u2.y*Math.sin(w);
		Epp.z = Hpp.z/2 + h*u1.z*Math.cos(w) + h*u2.z*Math.sin(w);

		var Ep = {};
		Ep.x = Epp.x*Math.cos(th) + Epp.z*Math.sin(th);
		Ep.y = Epp.y;
		Ep.z = Epp.x*Math.sin(th) + Epp.z*Math.cos(th);	

		var E = {};
		E.x = Ep.x + S.x;
		E.y = Ep.y + S.y;
		E.z = Ep.z + S.z;


		return E;
	}	

	function getPropQuaternion (tossOrientation, rotationAxis, jugglerRotation, propRotation, hand) {

		T = new THREE.Vector3(tossOrientation.x,tossOrientation.y,tossOrientation.z);
		R = new THREE.Vector3(rotationAxis.x,rotationAxis.y,rotationAxis.z);
		C = new THREE.Vector3(0,-1,0);

		if (hand == LEFT) {
			T.x *= -1;	
			R.x *= -1;
		}

		// rotate by juggler's rotation
		var Q1 = new THREE.Quaternion();
		Q1.setFromAxisAngle(new THREE.Vector3(0,-1,0), jugglerRotation);
		T.applyQuaternion(Q1);
		
		// get rotation to tossOrientation
		var Q2 = new THREE.Quaternion();
		Q2.setFromAxisAngle(new THREE.Vector3(T.z,0,-T.x), Math.acos(T.y));		
		
		// rotate rotationAxis according to tossOrientation
		var RQ = new THREE.Quaternion();
		RQ.setFromAxisAngle(new THREE.Vector3(0,1,0), Math.acos(-R.z*T.x + R.x*T.z));
		R.applyQuaternion(RQ);
		
		// rotate according to prop rotation
		var Q3 = new THREE.Quaternion();
		Q3.setFromAxisAngle(R, propRotation);
		
		// return composite rotation
		var q = new THREE.Quaternion();
		q = (new THREE.Quaternion()).multiplyQuaternions(Q2,Q3);
		
		return q;

	}

	function checkForCollision() {
		var r1, r2;
		// iterate over each props positions array
		for (var i = 0; i < siteswap.propPositions.length; i++) {			
			r1 = siteswap.props[i].radius;
			// iterate over all positions
			for (var j = 0; j < siteswap.propPositions[i].length; j++) {				
				// check position against all other props at that time
				for (var k = i+1; k < siteswap.propPositions.length; k++) {
					r2 = siteswap.props[i].radius;
					if (
						Math.sqrt(
							Math.pow(siteswap.propPositions[i][j].x-siteswap.propPositions[k][j].x,2)+
							Math.pow(siteswap.propPositions[i][j].y-siteswap.propPositions[k][j].y,2)+
							Math.pow(siteswap.propPositions[i][j].z-siteswap.propPositions[k][j].z,2)
						) <= (r1+r2)
					) {
						return true;
					}
				}
			}
		}
		return false;
	}

}

module.exports.CreateSiteswap = CreateSiteswap;
},{"./Bezier.js":1,"./BounceGA.js":2,"./util.js":7}],4:[function(require,module,exports){
module.exports.SiteswapAnimator = function(containerId, options) {
	
	var 
		container,
		width,
		height,
		camera, 
		scene, 
		renderer,
		/* camera starting point */
		camTheta = Math.PI, 
		camPhi = .7, 
		camRadius = 5,
		/* helpers for mouse interaction */
		isMouseDown = false, 
		onMouseDownTheta, 
		onMouseDownPhi, 
		onMouseDownPosition,
		cameraMode = {mode: 'sky'},
		propMeshes = [],
		jugglerMeshes = [],
		surfaceMeshes = [],
		propPathLines = [],
		jugglerHandVertices,
		jugglerElbowVertices,
		highestPoint,
		startTime,
		siteswap,
		renderMode = '3D',
		context,
		randomColors = ['red','blue','green','black','yellow','purple'],
		drawHands = false;

	this.displayPropPaths = options.displayPropPaths === true ? true : false;

	container = $('#' + containerId);

	width = container.width();
	height = container.height();

	if (renderMode == '2D') {

		container.append('<canvas id="siteswapAnimatorCanvas"></canvas>');
		canvas = container.find('canvas')[0];
		canvas.height = height;
		canvas.width = width;
		context = canvas.getContext('2d');

		context.clearRect(0,0,width,height);

	} else if (renderMode == '3D') {
		
		camera = new THREE.PerspectiveCamera( 90, width / height, .05, 100 );
		updateCamera();

		scene = new THREE.Scene();		
		
		/* create the renderer and add it to the canvas container */
		/* if browser is mobile, render using canvas */
		if( !window.WebGLRenderingContext ) {
			renderer = new THREE.CanvasRenderer();	
		} else {
			renderer = new THREE.WebGLRenderer( {antialias: true, preserveDrawingBuffer: true} );
		}
		
		renderer.setSize( width, height );

		container.empty();
		container.append(renderer.domElement);

		//add the event listeners for mouse interaction
		renderer.domElement.addEventListener( 'mousemove', onDocumentMouseMove, false );
		renderer.domElement.addEventListener( 'touchmove', onDocumentTouchMove, false );
		renderer.domElement.addEventListener( 'mousedown', onDocumentMouseDown, false );
		renderer.domElement.addEventListener( 'touchstart', onDocumentTouchStart, false );
		renderer.domElement.addEventListener( 'mouseup', onDocumentMouseUp, false );
		renderer.domElement.addEventListener( 'touchend', onDocumentMouseUp, false );
		renderer.domElement.addEventListener( 'mousewheel', onDocumentMouseWheel, false );

		onMouseDownPosition = new THREE.Vector2();

		renderer.setClearColor( 0xffffff, 1);

		renderer.render(scene,camera);

		this.renderer = renderer; // expose renderer		
		this.oneCycleComplete = false;
		this.paused = false;
		this.animationSpeed = .6;

	}

	this.resize = function(w,h) {
		width = w;
		height = h;
		camera = new THREE.PerspectiveCamera( 75, width / height, .05, 100 );
		updateCamera();
		renderer.setSize(width, height);
	}

	this.init = function(s,options) {

		this.paused = false;
		startTime = 0;
		siteswap = s;

		if (siteswap.errorMessage) {
			
			$('#errorMessage').html("ERROR: " + siteswap.errorMessage);
			$('#errorMessage').show();

		} else {

			/* show warnings for doing passing/bouncing with rings/clubs */
			if (siteswap.pass) {
				for (var i = 0; i < siteswap.props.length; i++) {
					if (siteswap.props[i].type == 'club' || siteswap.props[i].type == 'club') {
						$('#errorMessage').html("WARNING: Passing patterns with clubs/rings may look weird. Still working out kinks with prop orientation.");
						$('#errorMessage').show();
						break;
					}
				}
			}

			/* find highest point in the pattern */
			highestPoint = 0;
			for (var i = 0; i < siteswap.propPositions.length; i++) {
				for (var j = 0; j < siteswap.propPositions[i].length; j++) {
					if (siteswap.propPositions[i][j].y > highestPoint) {
						highestPoint = siteswap.propPositions[i][j].y;
					}
				}
			}
			camRadius = highestPoint+.5;

			/* clear out all meshes from scene */
			for( var i = scene.children.length - 1; i >= 0; i--) { scene.remove(scene.children[i]); }

			/* lights */
			var ceilingLight = new THREE.PointLight( 0xffffff );
			ceilingLight.position.set(0,20,0);
			scene.add( ceilingLight );
			var floorLight = new THREE.PointLight( 0xffffff );
			floorLight.position.set(0,0,-2);
			scene.add( floorLight );

			propMeshes = [];
			jugglerMeshes = [];
			jugglerHandVertices = [];
			jugglerElbowVertices = [];
			propPathLines = [];

			drawHands = options.drawHands ? options.drawHands : 0;

			drawSurfaces();

			drawJugglers();

			/* create each prop and add to empty propMeshes array */
			for (var i = 0; i < siteswap.numProps; i++) {

				var geometry;

				if (siteswap.props[i].type == "ball") {
					geometry = new THREE.SphereGeometry( siteswap.props[i].radius, 20 );
				}
				else if (siteswap.props[i].type == "club") {
					geometry = new THREE.CylinderGeometry( .008, .02, .02, 7, 5 );
					geometry.vertices.map(function(v) { v.y += .01; });
					var clubHandle = new THREE.CylinderGeometry( .015, .008, .18, 7, 5 );
					clubHandle.vertices.map(function(v) { v.y += .11; });
					var clubBody1 = new THREE.CylinderGeometry( .04, .015, .18, 7, 5 );
					clubBody1.vertices.map(function(v) { v.y += .29});
					var clubBody2 = new THREE.CylinderGeometry( .02, .04, .11, 7, 5 );
					clubBody2.vertices.map(function(v) { v.y += .43});
					THREE.GeometryUtils.merge(geometry,clubHandle);
					THREE.GeometryUtils.merge(geometry,clubBody1);
					THREE.GeometryUtils.merge(geometry,clubBody2);
					// move entire club down to correct center of gravity
					geometry.vertices.map(function(v) { v.y -= .2});

				}
				else if (siteswap.props[i].type == "ring") {
					// ring meshes
					var points = [];
					points.push( new THREE.Vector3( .14, 0, .01 ) );
					points.push( new THREE.Vector3( .18, 0, .01 ) );
					points.push( new THREE.Vector3( .18, 0, -.01 ) );
					points.push( new THREE.Vector3( .14, 0, -.01 ) );
					points.push( new THREE.Vector3( .14, 0, .01 ) );
					geometry = new THREE.LatheGeometry( points );

				}

				if (options.motionBlur) {
					var numTails = 2;
				} else {
					var numTails = 0;
				}
				
				var tmpPropMeshes = [];
				
				var propColor = (siteswap.props[i].color == "random" ? randomColors[Math.floor(Math.random()*randomColors.length)] : siteswap.props[i].color);

				for (var j = 0; j <= numTails; j++) {

					var material;					

					if (j == 0) {
						material = new THREE.MeshLambertMaterial( { color: propColor } );
						//material = new THREE.MeshBasicMaterial( { color: propColor, wireframe: true } );
					} else {
						material = new THREE.MeshLambertMaterial( { color: propColor, transparent: true, opacity: 1-1/(numTails+1)*j } );
					}
					var mesh = new THREE.Mesh( geometry, material );			

					scene.add( mesh );

					tmpPropMeshes.push(mesh); 

				}

				propMeshes.push( tmpPropMeshes );
				
			}

			if (this.displayPropPaths) {
				buildPropPaths();			
			}

		}

	}

	this.animate = function() {		

		if (startTime === 0) {
			startTime = (new Date()).getTime();
		}

		/* find time in the pattern and translate that to a discrete step in the prop position arrays */
		var timeElapsed = ((new Date()).getTime() - startTime)*this.animationSpeed;
		if (timeElapsed > (siteswap.states.length*siteswap.beatDuration*1000)) {
			this.oneCycleComplete = true;
		}
		var t = timeElapsed % (siteswap.states.length*siteswap.beatDuration*1000); // need to *1000 b/c timeElapsed is in ms
		
		this.render(t);

		if (!this.paused) {
			var self = this;
			requestAnimationFrame(function() { self.animate(); });
		}

	}

	this.render = function(t) {

		var step = Math.floor(t/(siteswap.states.length*siteswap.beatDuration*1000)*siteswap.numSteps);

		/* update prop mesh positions and rotations */
		for (var i = 0; i < propMeshes.length; i++) {
			for (var j = 0; j < propMeshes[i].length; j++) {

				var stepIx = step-j*Math.floor(siteswap.numStepsPerBeat/8); // the 10 here is the tail length factor
				if (stepIx < 0) {
					stepIx += siteswap.numSteps; 
				}

				propMeshes[i][j].position.x = siteswap.propPositions[i][stepIx].x;
				propMeshes[i][j].position.y = siteswap.propPositions[i][stepIx].y;
				propMeshes[i][j].position.z = siteswap.propPositions[i][stepIx].z;

				/* apply current rotation */				
				propMeshes[i][j].quaternion.set(1,0,0,0);
				
				// rotate rings so they are in correct position by default
				if (siteswap.props[i].type == 'ring') {
					var rotateRing = new THREE.Quaternion();
					rotateRing.setFromAxisAngle(new THREE.Vector3(0,1,0), Math.PI/2);
					propMeshes[i][j].quaternion.multiply(rotateRing);
				}

				var q = siteswap.propRotations[i][stepIx];
				propMeshes[i][j].quaternion.multiplyQuaternions(q, propMeshes[i][j].quaternion);

			}
		}

		/* update juggler hand positions */
		updateHandAndElbowPositions(step);

		updateCamera();

		/* mark geometry vertices as needs update */
		for (var i = 0; i < jugglerMeshes.length; i++) {
			for (var j = 0; j < jugglerMeshes[i].children.length; j++) {
				jugglerMeshes[i].children[j].geometry.verticesNeedUpdate = true;
			} 
		}

		try {
			renderer.render(scene, camera);
		} catch(e) {
			console.log('Error rendering');
		}

	}

	function drawSurfaces() {

		siteswap.surfaces.map(function(a) {
			var surface = {
				position: new THREE.Vector3(a.position.x,a.position.y,a.position.z),
				normal: new THREE.Vector3(a.normal.x,a.normal.y,a.normal.z),
				axis1: new THREE.Vector3(a.axis1.x,a.axis1.y,a.axis1.z),
				axis2: new THREE.Vector3(a.axis2.x,a.axis2.y,a.axis2.z),
				scale: a.scale
			}
			var surfaceGeom = new THREE.Geometry();
			surfaceGeom.vertices.push( (new THREE.Vector3()).copy(surface.position).add((new THREE.Vector3()).add(surface.axis1).add(surface.axis2)) );
			surfaceGeom.vertices.push( (new THREE.Vector3()).copy(surface.position).add((new THREE.Vector3()).add(surface.axis1).negate().add(surface.axis2)) );
			surfaceGeom.vertices.push( (new THREE.Vector3()).copy(surface.position).add((new THREE.Vector3()).add(surface.axis1).add(surface.axis2).negate()) );
			surfaceGeom.vertices.push( (new THREE.Vector3()).copy(surface.position).add((new THREE.Vector3()).add(surface.axis2).negate().add(surface.axis1)) );

			surfaceGeom.faces.push( new THREE.Face3( 0, 1, 2 ) );
			surfaceGeom.faces.push( new THREE.Face3( 2, 0, 3 ) );

			var surfaceMesh = new THREE.Mesh(surfaceGeom, new THREE.MeshBasicMaterial( { color: 'grey', side: THREE.DoubleSide } ));
			surfaceMeshes.push(surfaceMesh);
			scene.add(surfaceMesh);
		});
	}

	function toVector3(v) {
		return new THREE.Vector3(v.x,v.y,v.z);
	}

	function drawJugglers() {

		/* create each juggler and add to empty jugglerMeshes array */
		for (var i = 0; i < siteswap.numJugglers; i++) {

			function transformVector(x,y,z) {
				return new THREE.Vector3(
					siteswap.jugglers[i].position.x+(Math.cos(siteswap.jugglers[i].rotation)*x-Math.sin(siteswap.jugglers[i].rotation)*z),
					y,
					siteswap.jugglers[i].position.z+(Math.cos(siteswap.jugglers[i].rotation)*z+Math.sin(siteswap.jugglers[i].rotation)*x)
				);
			}

			jugglerHandVertices.push([[],[]]);
			jugglerElbowVertices.push([[],[]]);
			
			/* create juggler mesh at 0,0,0 */

			var jugglerMesh = new THREE.Object3D();

			var jugglerLegsG = new THREE.Geometry();
			jugglerLegsG.vertices.push(transformVector(-.125,0,0));
			jugglerLegsG.vertices.push(transformVector(0,.8,0));
			jugglerLegsG.vertices.push(transformVector(.125,0,0));
			var jugglerLegs = new THREE.Line(jugglerLegsG, new THREE.LineBasicMaterial({color: 'black'}));

			var jugglerTorsoG = new THREE.Geometry();
			jugglerTorsoG.vertices.push(transformVector(0,.8,0));
			jugglerTorsoG.vertices.push(transformVector(0,1.5,0));
			var jugglerTorso = new THREE.Line(jugglerTorsoG, new THREE.LineBasicMaterial({color: 'black'}));

			var jugglerShouldersG = new THREE.Geometry();
			jugglerShouldersG.vertices.push(transformVector(-.225,1.425,0));
			jugglerShouldersG.vertices.push(transformVector(.225,1.425,0));
			var jugglerShoulders = new THREE.Line(jugglerShouldersG, new THREE.LineBasicMaterial({color: 'black'}));	

			for (var h = -1; h <= 1; h+=2) {
				
				var hix = (h == -1 ? 0 : 1);
				
				armG = new THREE.Geometry();
				armG.vertices.push(transformVector(h*.225,1.425,0));
				jugglerElbowVertices[i][hix] = transformVector(h*.225,1.0125,0);
				
				armG.vertices.push(jugglerElbowVertices[i][hix]);	
				jugglerHandVertices[i][hix].push(transformVector(h*.225,1.0125,-.4125));
				armG.vertices.push(jugglerHandVertices[i][hix][0]);
				armG.dynamic = true;
				var arm = new THREE.Line(armG, new THREE.LineBasicMaterial({color: 'black'}));

				jugglerMesh.add( arm );

				if (drawHands) {

					// hand, don't need to specify vertices now because they'll be set later
					var handG = new THREE.Geometry();
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					handG.vertices.push(jugglerHandVertices[i][hix][1]);
					handG.vertices.push(jugglerHandVertices[i][hix][2]);
					handG.vertices.push(jugglerHandVertices[i][hix][3]);
					handG.vertices.push(jugglerHandVertices[i][hix][4]);
					handG.faces.push( new THREE.Face3( 0, 1, 2 ) );
					handG.faces.push( new THREE.Face3( 2, 0, 3 ) );
					handG.dynamic = true;
					var handMesh = new THREE.Mesh(handG, new THREE.MeshBasicMaterial( { color: 'black', side: THREE.DoubleSide }));

					jugglerMesh.add( handMesh );

					var fingerG = new THREE.Geometry();
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					fingerG.vertices.push(jugglerHandVertices[i][hix][5]);
					fingerG.vertices.push(jugglerHandVertices[i][hix][6]);
					fingerG.vertices.push(jugglerHandVertices[i][hix][7]);
					fingerG.dynamic = true;
					var finger = new THREE.Line(fingerG, new THREE.LineBasicMaterial({color: 'black'}));

					jugglerMesh.add( finger );

					fingerG = new THREE.Geometry();
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					fingerG.vertices.push(jugglerHandVertices[i][hix][8]);
					fingerG.vertices.push(jugglerHandVertices[i][hix][9]);
					fingerG.vertices.push(jugglerHandVertices[i][hix][10]);
					fingerG.dynamic = true;
					finger = new THREE.Line(fingerG, new THREE.LineBasicMaterial({color: 'black'}));

					jugglerMesh.add( finger );

					fingerG = new THREE.Geometry();
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					fingerG.vertices.push(jugglerHandVertices[i][hix][11]);
					fingerG.vertices.push(jugglerHandVertices[i][hix][12]);
					fingerG.vertices.push(jugglerHandVertices[i][hix][13]);
					fingerG.dynamic = true;
					finger = new THREE.Line(fingerG, new THREE.LineBasicMaterial({color: 'black'}));

					jugglerMesh.add( finger );

					fingerG = new THREE.Geometry();
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					jugglerHandVertices[i][hix].push(new THREE.Vector3());
					fingerG.vertices.push(jugglerHandVertices[i][hix][14]);
					fingerG.vertices.push(jugglerHandVertices[i][hix][15]);
					fingerG.vertices.push(jugglerHandVertices[i][hix][16]);
					fingerG.dynamic = true;
					finger = new THREE.Line(fingerG, new THREE.LineBasicMaterial({color: 'black'}));

					jugglerMesh.add( finger );

				}


			}

			var jugglerHead = new THREE.Mesh( new THREE.SphereGeometry( .1125, 20 ), new THREE.MeshPhongMaterial( { color: 'black' } ) );
			jugglerHead.position = new THREE.Vector3(siteswap.jugglers[i].position.x+Math.cos(siteswap.jugglers[i].rotation)*0,1.6125,siteswap.jugglers[i].position.z+Math.sin(siteswap.jugglers[i].rotation)*0);
			
			jugglerMesh.add( jugglerLegs );
			jugglerMesh.add( jugglerTorso );
			jugglerMesh.add( jugglerShoulders );
			jugglerMesh.add( jugglerHead );

			scene.add(jugglerMesh);
			jugglerMeshes.push(jugglerMesh);

		}
	}

	function buildPropPaths() {

		for (var i = 0; i < siteswap.propPositions.length; i++) {
			var propPathGeom = new THREE.Geometry();
			for (var j = 0; j < siteswap.propPositions[i].length; j++) {
				var propPosition = siteswap.propPositions[i][j];
				var eps = .001;
				propPathGeom.vertices.push(new THREE.Vector3(propPosition.x+(Math.random()-.5)*eps,propPosition.y+(Math.random()-.5)*eps,propPosition.z+(Math.random()-.5)*eps));
			}
			var propPathLine = new THREE.Line(propPathGeom, new THREE.LineBasicMaterial({color: siteswap.props[i].color}));
			propPathLines.push(propPathLine);
			scene.add(propPathLine);
		}

	}

	this.hidePropPaths = function() {
		propPathLines.map(function(a) { a.visible = false; });
	}

	this.showPropPaths = function() {
		if (propPathLines.length == 0) {
			buildPropPaths();
		} else {
			propPathLines.map(function(a) { a.visible = true; });
		}		
	}

	function updateHandAndElbowPositions(step) {
				
		for (var i = 0; i < jugglerHandVertices.length; i++) {
			for (var j = 0; j < 2; j++) {
				
				if (drawHands) {
					var handSize = .03;
					var zOffset = .03;
					var handVerticesDiff = [
						new THREE.Vector3(0,0,handSize+zOffset),
						new THREE.Vector3(-handSize,0,-handSize+zOffset),
						new THREE.Vector3(handSize,0,-handSize+zOffset),
						new THREE.Vector3(handSize,0,handSize+zOffset),
						new THREE.Vector3(-handSize,0,handSize+zOffset),
						new THREE.Vector3(0,0,-handSize+zOffset),
						new THREE.Vector3(0,.01,-2*handSize+zOffset),
						new THREE.Vector3(0,.05,-3*handSize+zOffset),
						new THREE.Vector3(handSize,0,-handSize+zOffset),
						new THREE.Vector3(handSize,.01,-2*handSize+zOffset),
						new THREE.Vector3(handSize,.05,-3*handSize+zOffset),
						new THREE.Vector3(-handSize,0,-handSize+zOffset),
						new THREE.Vector3(-handSize,.01,-2*handSize+zOffset),
						new THREE.Vector3(-handSize,.05,-3*handSize+zOffset),
						new THREE.Vector3((j == 0 ? -1 : 1)*handSize,0,handSize+zOffset),
						new THREE.Vector3((j == 0 ? -1 : 1)*2*handSize,0,handSize+zOffset-.02),
						new THREE.Vector3((j == 0 ? -1 : 1)*2*handSize,.02,handSize+zOffset-.04)
					];
					var angle = siteswap.jugglerHandPositions[i][j][step].angle;
					var propRadius = siteswap.props[0].radius;
					var jugglerWristPosition = toVector3(siteswap.jugglerHandPositions[i][j][step]).add(new THREE.Vector3(propRadius*Math.sin(angle),-propRadius*Math.cos(angle),0));
					
					var handAngle = -Math.atan2(siteswap.jugglerElbowPositions[i][j][step].x-siteswap.jugglerHandPositions[i][j][step].x,siteswap.jugglerElbowPositions[i][j][step].z-siteswap.jugglerHandPositions[i][j][step].z);

					for (var k = 0; k < handVerticesDiff.length; k++) {
						var newX = handVerticesDiff[k].x*Math.cos(angle) - handVerticesDiff[k].y*Math.sin(angle);
						var newY = handVerticesDiff[k].y*Math.cos(angle) + handVerticesDiff[k].x*Math.sin(angle);
						newX = newX*Math.cos(handAngle) - handVerticesDiff[k].z*Math.sin(handAngle);
						var newZ = handVerticesDiff[k].z*Math.cos(handAngle) + newX*Math.sin(handAngle);
						handVerticesDiff[k].x = newX;
						handVerticesDiff[k].y = newY;
						handVerticesDiff[k].z = newZ;
						jugglerHandVertices[i][j][k].copy((new THREE.Vector3()).copy(jugglerWristPosition).add(handVerticesDiff[k]));
					}
				} else {
					jugglerHandVertices[i][j][0].copy(toVector3(siteswap.jugglerHandPositions[i][j][step])); 
				}

			}

			jugglerElbowVertices[i][0].x = siteswap.jugglerElbowPositions[i][0][step].x;
			jugglerElbowVertices[i][0].y = siteswap.jugglerElbowPositions[i][0][step].y;
			jugglerElbowVertices[i][0].z = siteswap.jugglerElbowPositions[i][0][step].z;
			jugglerElbowVertices[i][1].x = siteswap.jugglerElbowPositions[i][1][step].x;
			jugglerElbowVertices[i][1].y = siteswap.jugglerElbowPositions[i][1][step].y;
			jugglerElbowVertices[i][1].z = siteswap.jugglerElbowPositions[i][1][step].z;

		}
	}

	function updateCamera() {
		if (cameraMode.mode == 'sky') {
			camera.position.x = camRadius * Math.sin( camTheta ) * Math.cos( camPhi );
			camera.position.y = camRadius * Math.sin( camPhi );
			camera.position.z = camRadius * Math.cos( camTheta ) * Math.cos( camPhi );
			var lookAt = new THREE.Vector3(0,0,0);
			if (siteswap !== undefined) {
				for (var i = 0; i < siteswap.jugglers.length; i++) {
					lookAt.x += siteswap.jugglers[i].position.x;
					lookAt.z += siteswap.jugglers[i].position.z;
				}
				lookAt.x /= siteswap.jugglers.length;
				lookAt.z /= siteswap.jugglers.length;
				lookAt.y = highestPoint/2;
			}			
			camera.lookAt(lookAt);
		} else if (cameraMode.mode == 'juggler') {
			/* need to update x and y to reflect the position of the juggler you are possessing */
			camera.position.x = 0;
			camera.position.y = 1.6125;
			camera.position.z = 0;
			//camera.lookAt(new THREE.Vector3(Math.sin(camTheta),3,Math.cos(camTheta)));
			camera.lookAt(new THREE.Vector3(Math.sin(camTheta)*Math.cos(camPhi),1.6125-Math.sin(camPhi),Math.cos(camTheta)*Math.cos(camPhi)));
		} else if (cameraMode.mode == 'custom') {
			var cameraPosition = new THREE.Vector3(cameraMode.x,cameraMode.y,cameraMode.z);
			camera.position = cameraPosition;
			camera.lookAt(new THREE.Vector3(cameraPosition.x-Math.sin(camTheta)*Math.cos(camPhi),cameraPosition.y-Math.sin(camPhi),cameraPosition.z-Math.cos(camTheta)*Math.cos(camPhi)));
		}
	}

	this.zoomIn = function() { camRadius-=.1; }

	this.zoomOut = function() { camRadius+=.1; }

	/* got the camera rotation code from: http://www.mrdoob.com/projects/voxels/#A/ */
	function onDocumentMouseDown( event ) {
		isMouseDown = true;
		onMouseDownTheta = camTheta;
		onMouseDownPhi = camPhi;
		onMouseDownPosition.x = event.clientX;
		onMouseDownPosition.y = event.clientY;
	}

	function onDocumentTouchStart( event ) {
		isMouseDown = true;
		onMouseDownTheta = camTheta;
		onMouseDownPhi = camPhi;
		onMouseDownPosition.x = event.changedTouches[0].clientX;
		onMouseDownPosition.y = event.changedTouches[0].clientY;
	}

	function onDocumentMouseMove( event ) {
		event.preventDefault();
		if ( isMouseDown ) {
			camTheta = - ( ( event.clientX - onMouseDownPosition.x ) * 0.01 ) + onMouseDownTheta;
			
			var dy = event.clientY - onMouseDownPosition.y;
			
			var newCamPhi = ( ( dy ) * 0.01 ) + onMouseDownPhi;

			if (newCamPhi < Math.PI/2 && newCamPhi > -Math.PI/2) {
				camPhi = newCamPhi;
			}
		}

		updateCamera();
		renderer.render(scene, camera);
	}

	function onDocumentTouchMove( event ) {
		event.preventDefault();
		if ( isMouseDown ) {
			camTheta = - ( ( event.changedTouches[0].clientX - onMouseDownPosition.x ) * 0.01 ) + onMouseDownTheta;
			
			var dy = event.changedTouches[0].clientY - onMouseDownPosition.y;
			
			var newCamPhi = ( ( dy ) * 0.01 ) + onMouseDownPhi;

			if (newCamPhi < Math.PI/2 && newCamPhi > -Math.PI/2) {
				camPhi = newCamPhi;
			}
		}

		updateCamera();
		renderer.render(scene, camera);
	}

	function onDocumentMouseUp( event ) {
		event.preventDefault();
		isMouseDown = false;
	}

	function onDocumentMouseWheel( event ) { camRadius -= event.wheelDeltaY*.002; }

	this.updateAnimationSpeed = function(speed) {
		this.animationSpeed = speed;
	}

	this.updateCameraMode = function(mode) {
		cameraMode = mode;
	}

}
},{}],5:[function(require,module,exports){
// some helper functions 
if (!Array.prototype.last){
    Array.prototype.last = function(){
        return this[this.length - 1];
    };
};

if (!Array.prototype.sum){
    Array.prototype.sum = function(){
        return this.reduce(function(a,b) { return a+b; });
    };
};

function factorial(a) { if (a == 2) { return a; } else if (a == 0) { return 1; } else { return a*factorial(a-1); } } 

module.exports.siteswapGraph = function(config, outputs) {
	
	// apply default configs 
	config.minPeriod = (config.minPeriod === undefined ? [] : config.minPeriod);
	config.async = (config.async === undefined ? [] : config.async);
	config.callbacks = (config.callbacks === undefined ? [] : config.callbacks);
	config.exclude = (config.exclude === undefined ? [] : config.exclude);
	config.includeExcited = (config.includeExcited === undefined ? true : config.includeExcited);
	config.maxSearches = (config.maxSearches === undefined ? 99999999 : config.maxSearches);
	config.maxSiteswaps = (config.maxSiteswaps === undefined ? 1000 : config.maxSiteswaps);
	config.includeMultiplex = (config.includeMultiplex === undefined ? false : config.includeMultiplex);
	config.sync = (config.sync === undefined ? false : config.sync);
	
	// init outputs object, this will contain the graph and will be updated as we go for async access 
	outputs = (outputs === undefined ? {} : outputs);

	var nodes = [];			// each node is a state (ie. prop landing schedule)
	var edges = [];			// each edge specifies a source/target node and a transition value (ie. the siteswap)
	var siteswaps = [];		// each siteswap is an array of edges that composes a unique cycle
	var formattedSiteswaps = [];

	outputs.graph = {
		nodes: nodes,
		edges: edges,
		siteswaps: siteswaps,
		formattedSiteswaps: formattedSiteswaps
	};


	// build nodes helper function, this also kicks off building the edges once all nodes have been created 
	function buildNodes(node,nodeOptions,last) {

		// if the node we're constructing has reached the expected length
		if (node.length == config.maxPeriod) {
			if (node.reduce(function(a,b) { return a+b; }) == config.numProps) {
				nodes.push({value: node, edges: []});
			}			
			
			if (config.callbacks.updateGraphProgress) {				
				//var progress = nodes.length/expectedNumNodes;
				var progress = 1; // todo
				config.callbacks.updateGraphProgress(progress);
			}	
			// if we've created all the nodes kick off the function to build the edges
			if (last) {
				buildEdges();
			}

		// if the node is not the expected length then we need to keep building
		} else {

			nodeOptions.map(function(nodeOption,ix,nodeOptions) {		

				var newNode = node.slice();
				newNode.push(nodeOption);

				var propDiff = newNode.reduce(function(a,b) { return a+b; }) - config.numProps;
				var newNodeOptions = [];
				
				// always have 1 first so the first node is the ground node
				if (propDiff == 0) {
					newNodeOptions = [0];
				} else if (propDiff == 1) {
					newNodeOptions = [1,0];
				} else {
					if (config.includeMultiplex) {
						newNodeOptions = [1,2,0];
					} else {
						newNodeOptions = [1,0];
					}
					
				}

				var last = this.last && (ix == nodeOptions.length-1);
				if (config.async) {
					setTimeout(buildNodes,0,newNode,newNodeOptions,last);
				} else {
					buildNodes(newNode,newNodeOptions,last);
				}
			
			},{last:last}); // no idea why i have to do this, for some reason the variable last isn't available within the scope of the map call. but the variable node is. wtf?

		}
	}

	// helper functions to get edges between 2 nodes
	function getEdgesBetween2Nodes(node1,node2) {

		var edges = [];

		var nextUp = node1[0];
		var newNode = node1.slice(1,node1.length);
		newNode.push(0);
		
		var multiplex = false;
		if (nextUp == 0) {
			edges.push('0');
		} else if (nextUp > 1) {
			multiplex = true;
			edges.push('[');
		} else {
			edges.push('');
		}
		for (var i = 0; i < newNode.length; i++) {
			var tossValue = i+1;
			if (tossValue > 9) {
				tossValue = String.fromCharCode(87+tossValue);
			}
			if(newNode[i] != node2[i]) {
				if (nextUp >= (node2[i] - newNode[i])) {
					edges[0] += tossValue;
					nextUp--;
					if (nextUp == 1 && (node2[i] - newNode[i]) == 2) {
						edges[0] += tossValue;
						nextUp--;
					}	
				} else {
					return [];
				}
			}
		}
		if (multiplex) {
			edges[0] += ']';
		}	

		return edges;
	}

	// recursively build edges between all nodes, want to make this async
	function buildEdges() {

		// compare all nodes
		for (var i = 0; i < nodes.length; i++) {
			if (config.callbacks.updateGraphProgress) {
				var progress = (i+1)/nodes.length;
				config.callbacks.updateGraphProgress(progress);
			}
			for (var j = 0; j < nodes.length; j++) {				
				function addEdgesToGraph(i,j) {
					// get edges between 2 nodes
					var edgeValues = getEdgesBetween2Nodes(nodes[i].value, nodes[j].value);
					// add each edge
					for (var k = 0; k < edgeValues.length; k++) {
						nodes[i].edges.push(edges.push({source: i, target: j, value: edgeValues[k]})-1);
					}
				}
				setTimeout(addEdgesToGraph,0,i,j);			
			}
		}

		if (config.callbacks.graphDone) {
			setTimeout(config.callbacks.graphDone);
		}

		// kick off siteswap search
		if (config.async) {
			setTimeout(search,0,0,[]);
		} else {
			search(0,[]);
		}
	}
	
	// compare 2 siteswap patterns to check for equality
	function patternsMatch(p1,p2) {

		if (p1.length != p2.length) {
			return false;
		} else {			
			for (var i = 0; i <= p1.length; i++) {
				if (p1.toString() == p2.toString()) {
					return true;
				}
				p1.push(p1[0]);
				p1 = p1.slice(1);	
				if (config.sync) {
					p1.push(p1[0]);
					p1 = p1.slice(1);
				}		
			}
		}
	}

	var numSearches = 0;	

	function search(origNodeIx,history) {

		numSearches++;	
		if (numSearches <= config.maxSearches && siteswaps.length < config.maxSiteswaps) {

			var nextNodeIx = origNodeIx;
			var validSiteswap = false;
			var siteswapStartNode = undefined;

			if (history.length > 0 && history.length <= config.maxPeriod) {
				// check if valid siteswap, ie. the last edge returns us to the first node
				if (history.length >= config.minPeriod && (!config.sync || history.length % 2 == 0)) {
					if (edges[history.last()].target == origNodeIx) {
						validSiteswap = true;
						siteswapStartNode = 0;
					} else {
						// excited siteswaps would return us to any node within the search (assuming the first node is the ground node)
						if (config.includeExcited) {
							for (var i = 0; i < history.length; i++) {
								if (edges[history.last()].target == edges[history[i]].source && (!config.sync || i % 2 == 0)) {
									validSiteswap = true;									
									siteswapStartNode = i;
								}
							}
						}					
					}
				}
				// if the siteswap is valid check to see if it exists or not, then add it to the list
				if (validSiteswap) {
					var siteswap = history.map(function (a,ix) {
						if (!config.sync) {
							return edges[a].value; 
						} else {
							var syncEdgeValue = "";
							var asyncEdgeValue = edges[a].value;
							for (var i = 0; i < asyncEdgeValue.length; i++) {
								if (asyncEdgeValue[i] == "[" || asyncEdgeValue[i] == "]" || parseInt(asyncEdgeValue[i]) % 2 == 0 ) {
									syncEdgeValue += asyncEdgeValue[i];
								} else if (ix % 2 == 0) {
									syncEdgeValue += ((parseInt(asyncEdgeValue[i])-1)+"x");
								} else {
									syncEdgeValue += ((parseInt(asyncEdgeValue[i])+1)+"x");
								}
							}
							return syncEdgeValue;
						}
					}).slice(siteswapStartNode);
					var exists = false;
					for (var i = 0; i < siteswaps.length; i++) {
						var ssToMatch = siteswaps[i].map(function (a,ix) {
							if (!config.sync) {
								return edges[a].value; 
							} else {
								var syncEdgeValue = "";
								var asyncEdgeValue = edges[a].value;
								for (var i = 0; i < asyncEdgeValue.length; i++) {
									if (asyncEdgeValue[i] == "[" || asyncEdgeValue[i] == "]" || parseInt(asyncEdgeValue[i]) % 2 == 0 ) {
										syncEdgeValue += asyncEdgeValue[i];
									} else if (ix % 2 == 0) {
										syncEdgeValue += ((parseInt(asyncEdgeValue[i])-1)+"x");
									} else {
										syncEdgeValue += ((parseInt(asyncEdgeValue[i])+1)+"x");
									}
								}
								return syncEdgeValue;
							}
						});
						if (patternsMatch(ssToMatch,siteswap.slice())) {							
							exists = true;
							break;
						}
					}
					if (!exists) {
						// the siteswaps array will actually store the edge history which can be converted into a siteswap string
						var siteswapIx = siteswaps.push(history.slice(siteswapStartNode))-1;
						if (config.callbacks.siteswapFound) {
							for (var i = 0; i < siteswap.length; i++) {
								if (siteswap[i] > 9) {
									siteswap[i] = String.fromCharCode(87+parseInt(siteswap[i]));
								}
							}
							var formattedSiteswap = "";
							if (config.sync) {
								for (var i = 0; i < siteswap.length; i++) {
									if (i % 2 == 0) {
										formattedSiteswap += ("(" + siteswap[i] + ",");
									} else { 
										formattedSiteswap += (siteswap[i] + ")");
									}
								}
							} else {
								formattedSiteswap = siteswap.join('');
							}
							formattedSiteswaps.push(formattedSiteswap);
							config.callbacks.siteswapFound(formattedSiteswap,siteswapIx,(siteswapStartNode > 0));							
						}
					}
					validSiteswap = true;
				}
				nextNodeIx = edges[history.last()].target;
			}

			// if the siteswap was invalid or it was valid but was excited, and we're still below the maxperiod for a siteswap, keep searching
			if ((!validSiteswap || siteswapStartNode > 0) && history.length < config.maxPeriod) {

				// search each edge of this next node
				nodes[nextNodeIx].edges.map(function(edgeIx) {

					function runSearch() {

						// was previously checking if we already visited the next node, but not doing that anymore since we're checking for excited swaps
						var alreadyVisited = false;
						// for (var j = 0; j < history.length; j++) {
						// 	if (edges[history[j]].target == edges[edgeIx].target) {
						// 		alreadyVisited = true;
						// 		break;
						// 	}
						// }

						// check if searching this edge is going to match the exclusion pattern
						// TODO: need to fix this to search better
						var exclude = false;
						for (var j = 0; j < config.exclude.length; j++) {
							if (config.exclude[j] == edges[edgeIx].value) {																
								exclude = true;
								break;
							}
						}

						// if this is an odd numbered edge in the history and we're doing sync, this can't be a 1
						if (config.sync && history.length % 2 == 0 && edges[edgeIx].value.indexOf(1) > -1) {
							exclude = true;
						}

						// execute the search through the edge
						if (!alreadyVisited && !exclude) {
							var newHistory = history.slice();
							newHistory.push(edgeIx);
							search(origNodeIx,newHistory);
						}

					}

					if (config.async) {
						setTimeout(runSearch,0);
					} else {
						runSearch();
					}

				});

			}
		}		
	}

	// kick off the whole process
	// get the values that can contribute to a state 
	var nodeOptions = [];
	nodeOptions.push(1);
	nodeOptions.push(0);
	if (config.includeMultiplex) {
		nodeOptions.push(2);
	}

	if (config.async) {
		setTimeout(buildNodes,0,[],nodeOptions.slice(),true);
	} else {			
		buildNodes([],nodeOptions.slice(),true);
	}

	// only explicitly return the outputs on synchronous mode. 
	// async calls should be monitoring the ouputs param
	if(!config.async) {
		return outputs;
	}

}
},{}],6:[function(require,module,exports){
var SiteswapJS = require('./Siteswap.js');
var SiteswapAnimator = require('./SiteswapAnimator.js');
var SiteswapGraph = require('./SiteswapGraph.js');
var util = require('./util.js');

var twoWindow = false;

window.onload = function () {

	displayMenu('pattern');

	updateAdvancedInputsFromBasic();	

	window.animator = new SiteswapAnimator.SiteswapAnimator('animatorCanvasContainer', {displayPropPaths: false});

	window.onresize();

	buildExamples();
	refreshSavedSiteswapsList();

	bindInputs(applyInputDefaults(getInputsFromQueryString()));

	go();

}

function displayMenu(menu) {	
	$('.controlDiv').hide()
	$('#'+menu+'Menu').show();
	$('#nav a').removeClass('selected');
	$('#nav a').addClass('unselected');
	$('#nav #' +menu).addClass('selected');
	$('#nav #' +menu).removeClass('unselected');	
	window.onresize();
}

window.onresize = function () {
	var windowWidth = $(window).width()-5;
	var windowHeight = $(window).height()-10;
	var controlsWidth = 500;
	var animatorWidth = windowWidth-controlsWidth;
	var minAnimatorWidth = 250;
	var animatorHeight = windowHeight;

	if (animatorWidth > minAnimatorWidth) {
		$('#nav #animator').hide();
		$('#animatorCanvasContainer').appendTo($('body'));
		$('#animatorMenu').removeClass('controlDiv');		
		twoWindow = true;
	} else {
		if (controlsWidth > windowWidth) {
			controlsWidth = windowWidth;
		}
		$('#nav #animator').show();
		$('#animatorCanvasContainer').appendTo($('#controlsContainer #animatorMenu'));
		$('#animatorMenu').addClass('controlDiv');
		animatorWidth = windowWidth;
		controlsWidth = windowWidth;
		twoWindow = false;
		animatorHeight = windowHeight - $('#animatorCanvasContainer').offset().top;
	}

	$('#controlsContainer').height(windowHeight);
	$('#controlsContainer').width(controlsWidth);

	$('#animatorCanvasContainer').height(animatorHeight);
	$('#animatorCanvasContainer').width(animatorWidth);	

	// resize divs containing lists
	$('#generatedSiteswaps').height(windowHeight-$('#generatedSiteswaps').offset().top);
	$('#exampleSiteswaps').height(windowHeight-$('#exampleSiteswaps').offset().top);
	$('#savedSiteswaps').height(windowHeight-$('#exampleSiteswaps').offset().top);
	$('#patternMenu').height(windowHeight-$('#patternMenu').offset().top);

	if (animator.resize) {
		animator.resize(animatorWidth, windowHeight);
	}	
}

function updateAdvancedInputsFromBasic() {
	bindInputs(applyInputDefaults({
		siteswap: $('#siteswap').val(),
		props: [{type: $('#prop').val(), color: 'red', radius: .05, C: .97}],
		beatDuration: $('#beatDuration').val(),
		dwellPath: $('#dwellPath').val()
	}));
	updateAdvancedLabels();
}

function applyInputDefaults(inputs) {
	inputs.siteswap = inputs.siteswap === undefined ? "3" : inputs.siteswap;
	inputs.props = inputs.props === undefined ? [{type: "ball", color: "red", radius: ".05", C: .97}] : inputs.props;
	inputs.beatDuration = inputs.beatDuration === undefined ? .28 : inputs.beatDuration;
	inputs.dwellRatio = inputs.dwellRatio === undefined ? .8 : inputs.dwellRatio;
	inputs.dwellPath = inputs.dwellPath === undefined ? "(30)(10)" : inputs.dwellPath;
	inputs.matchVelocity = inputs.matchVelocity === undefined ? 0 : inputs.matchVelocity;
	inputs.dwellCatchScale = inputs.dwellCatchScale === undefined ? .06 : inputs.dwellCatchScale;
	inputs.dwellTossScale = inputs.dwellTossScale === undefined ? .06 : inputs.dwellTossScale;
	inputs.emptyTossScale = inputs.emptyTossScale === undefined ? .025 : inputs.emptyTossScale;
	inputs.emptyCatchScale = inputs.emptyCatchScale === undefined ? .025 : inputs.emptyCatchScale;
	inputs.armAngle = inputs.armAngle === undefined ? .1 : inputs.armAngle;
	inputs.jugglers = inputs.jugglers === undefined ? [{position: {x: 0, z: 0}, rotation: 0}] : inputs.jugglers;
	inputs.surfaces = inputs.surfaces === undefined ? [{position: {x: 0, y: 0, z:0}, normal: {x: 0, y:1, z:0}, scale: 2}] : inputs.surfaces;
	return inputs;
}

function bindInputs(inputs) {
	var inputsText = inputs.siteswap + "\n";
	inputsText += inputs.beatDuration + " " + inputs.dwellRatio + "\n";
	for (var i = 0; i < inputs.props.length; i++) {
		inputsText += inputs.props[i].type + " " + inputs.props[i].color + " " + inputs.props[i].radius + " " + inputs.props[i].C;
		if (i < inputs.props.length-1) {
			inputsText += " ";
		} else {
			inputsText += "\n";
		}
	}
	inputsText += inputs.dwellPath + "\n";
	inputsText += inputs.matchVelocity + " " + inputs.dwellCatchScale + " " + inputs.dwellTossScale + " " + inputs.emptyTossScale + " " + inputs.emptyCatchScale + " " + inputs.armAngle + "\n";
	for (var i = 0; i < inputs.jugglers.length; i++) {
		inputsText += inputs.jugglers[i].position.x + " " + inputs.jugglers[i].position.z + " " + inputs.jugglers[i].rotation;
		if (i < inputs.jugglers.length-1) {
			inputsText += " ";
		}
	}
	if (inputs.surfaces.length > 0) {
		inputsText += "\n"
	}
	for (var i = 0; i < inputs.surfaces.length; i++) {
		inputsText += inputs.surfaces[i].position.x + " " + inputs.surfaces[i].position.y + " " + inputs.surfaces[i].position.z + " " + inputs.surfaces[i].normal.x + " " + inputs.surfaces[i].normal.y + " " + inputs.surfaces[i].normal.z + " " + inputs.surfaces[i].scale;
		if (i < inputs.surfaces.length-1) {
			inputsText += "\n";
		}
	}
	$('#inputsAdvanced').val(inputsText);
	updateAdvancedLabels();
} 

function parseInputs(inputs) {
	var lines = inputs.split('\n');
	var siteswap = lines[0];
	var beatDuration = parseFloat(lines[1].split(' ')[0]);
	var dwellRatio = parseFloat(lines[1].split(' ')[1]);
	var propsLine = lines[2].split(' ');
	var props = [];
	for (var i = 3; i < propsLine.length; i+=4) {
		props.push({
			type: propsLine[i-3],
			color: propsLine[i-2], 
			radius: parseFloat(propsLine[i-1]), 
			C: parseFloat(propsLine[i])
		});
	}
	
	// this whole bit should probably be moved into the Siteswap class
	var customDwellPathInput = lines[3];
	var customDwellPathBeats = customDwellPathInput.split(').').map(function(a,ix,arr) { if (ix < arr.length-1) { return a+')'; } else { return a; } });
	var dwellPath = [];
	for (var i = 0; i < customDwellPathBeats.length; i++) {
		var customDwellPathArr = customDwellPathBeats[i].match(/\(-?\d+(\.\d+)?(,-?\d+(\.\d+)?)?(,-?\d+(\.\d+)?)?(,\{-?\d+(\.\d+)?,-?\d+(\.\d+)?,-?\d+(\.\d+)?,-?\d+(\.\d+)?\})?\)/g);
		if ( customDwellPathArr.reduce(function(a,b) { return a+b }).length == customDwellPathBeats[i].length ) {
			dwellPath.push(
				customDwellPathArr.map(function(a,ix) {   
					var xyz = a.match(/\(-?\d+(\.\d+)?(,-?\d+(\.\d+)?)?(,-?\d+(\.\d+)?)?/g)[0].match(/-?\d+(\.\d+)?/g);
					var rot = a.match(/\{-?\d+(\.\d+)?,-?\d+(\.\d+)?,-?\d+(\.\d+)?,-?\d+(\.\d+)?\}/g); 
					var xyzth;
					if (rot) {
						xyzth = rot[0].match(/-?\d+(\.\d+)?/g);
					}
					var rotation;
					if (xyzth) {
						rotation = {x:parseFloat(xyzth[0]),y:parseFloat(xyzth[1]),z:parseFloat(xyzth[2]),th:parseFloat(xyzth[3])};
					} else if (props[0].type == 'club') {
						rotation = {x:4,y:0,z:(ix == 0 ? -1 : 1),th:Math.PI/2+(ix == 0 ? .5 : -.7)};
					} else if (props[0].type == 'ring') {
						rotation = {x:0,y:1,z:0,th:Math.PI/2};
					} else {
						rotation = {x:1,y:0,z:0,th:0};
					}
					return {
						x: parseFloat(xyz[0])/100,
						y: xyz[1] ? parseFloat(xyz[1])/100 : 0,
						z: xyz[2] ? parseFloat(xyz[2])/100 : 0,
						rotation: rotation
					}
				})
			);
		} else {
			throw 'Invalid custom dwell path';
		}
	}

	var dwellPathConfigs = lines[4].split(' ');
	var matchVelocity = dwellPathConfigs[0] == 1 ? true : false;
	var dwellCatchScale = parseFloat(dwellPathConfigs[1]);
	var dwellTossScale = parseFloat(dwellPathConfigs[2]);
	var emptyTossScale = parseFloat(dwellPathConfigs[3]);
	var emptyCatchScale = parseFloat(dwellPathConfigs[4]);
	var armAngle = parseFloat(dwellPathConfigs[5]);

	var jugglerPositions = lines[5].split(' ');
	var jugglers = undefined;
	for (var i = 0; i < jugglerPositions.length; i+=3) {
		if (jugglers === undefined) {
			jugglers = [];
		}
		jugglers.push({
			position: {x: parseFloat(jugglerPositions[i]), z: parseFloat(jugglerPositions[i+1])},
			rotation: parseFloat(jugglerPositions[i+2])
		});
	}

	var surfaces= [];
	for (var i = 6; i < lines.length; i++) {
		var surfaceLine = lines[i].split(' ');
		surfaces.push({
			position: {
				x: parseFloat(surfaceLine[0]),
				y: parseFloat(surfaceLine[1]),
				z: parseFloat(surfaceLine[2]),
			},
			normal: {
				x: parseFloat(surfaceLine[3]),
				y: parseFloat(surfaceLine[4]),
				z: parseFloat(surfaceLine[5]),	
			},
			scale: parseFloat(surfaceLine[6])
		});
	}

	return {
		siteswap: siteswap,
		beatDuration: beatDuration,
		dwellRatio: dwellRatio,
		props: props,
		inputDwellPath: customDwellPathInput,
		dwellPath: dwellPath,
		matchVelocity: matchVelocity,
		dwellCatchScale: dwellCatchScale,
		dwellTossScale: dwellTossScale,
		emptyTossScale: emptyTossScale,
		emptyCatchScale: emptyCatchScale,
		armAngle: armAngle,
		surfaces: surfaces,
		jugglers: jugglers
	};
}

function go() {

	if (!twoWindow) {
		displayMenu('animator');
	}

	var inputs = parseInputs($('#inputsAdvanced').val());

	var saveURL = window.location.href.replace("#","");
	if (saveURL.indexOf("?") > -1) {
		saveURL = saveURL.substring(0,saveURL.indexOf("?"));
	}

	var saveQueryString = "?v=16&siteswap=" + encodeURIComponent(inputs.siteswap) + "&beatDuration=" + inputs.beatDuration + "&dwellPath=" + inputs.inputDwellPath; 

	$('#saveURL').text(saveURL + saveQueryString);
	$('#saveURL').attr("href",saveURL + saveQueryString);	

	window.siteswap = SiteswapJS.CreateSiteswap(inputs.siteswap, 
		{
			beatDuration: 		inputs.beatDuration,
			dwellRatio: 		inputs.dwellRatio,
			props: 				inputs.props,
			dwellPath: 			inputs.dwellPath,
			matchVelocity: 		inputs.matchVelocity,
			dwellCatchScale: 	inputs.dwellCatchScale,
			dwellTossScale: 	inputs.dwellTossScale,
			emptyTossScale: 	inputs.emptyTossScale,
			emptyCatchScale: 	inputs.emptyCatchScale,
			armAngle: 			inputs.armAngle,
			surfaces: 			inputs.surfaces,
			jugglers: 			inputs.jugglers
		});

	if (siteswap.errorMessage) {
		animator.paused = true;
		$('#errorMessage').show();
		$('#message').text(siteswap.errorMessage);
	} else {

		if (siteswap.collision) {
			$('#errorMessage').show();
			$('#message').text("This pattern has collisions.");
		} else {
			$('#errorMessage').hide();
		}		

		var drawHands = false;
		if (siteswap.props[0].type == 'ball' && siteswap.numJugglers == 1) {
			drawHands = true;
		}

		animator.init(siteswap, 
			{
				drawHands: $('#drawHands')[0].checked
				//, motionBlur: true
			}
		);
		animator.animate();

	}

}

function zoomIn() { animator.zoomIn(); }

function zoomOut() { animator.zoomOut(); }

function updateAnimationSpeed() {
	var animationSpeed = parseFloat($('#animationSpeed').val());
	animator.updateAnimationSpeed(animationSpeed);
}

function updateCameraMode() {
	var mode = $('#cameraMode').val();
	cameraMode = {mode: mode};
	if (mode == "custom") {
		var cameraCustomPosition = $('#cameraCustomPosition').val().split(",");
		cameraMode.x = parseFloat(cameraCustomPosition[0]);
		cameraMode.y = parseFloat(cameraCustomPosition[1]);
		cameraMode.z = parseFloat(cameraCustomPosition[2]);
	}
	animator.updateCameraMode(cameraMode);
}

function updateDisplayPropPaths() {
	animator.displayPropPaths = !animator.displayPropPaths;
	if (animator.displayPropPaths) {
		animator.showPropPaths();
	} else {
		animator.hidePropPaths();
	}
}

function updateAdvancedLabels() {
	var inputs = parseInputs($('#inputsAdvanced').val());
	$('#lblSiteswap').text(inputs.siteswap);
	$('#lblBeatDuration').text(inputs.beatDuration);
	$('#lblDwellRatio').text(inputs.dwellRatio);
	var lblProps = "";
	for (var i = 0; i < inputs.props.length; i++) {
		var prop = inputs.props[i];
		lblProps += prop.color + ' ' + prop.type + ' ' + prop.radius + 'm ' + prop.C;
		if (i < inputs.props.length-1) {
			lblProps += ', ';
		}
	}
	$('#lblProps').text(lblProps);
	$('#lblDwellPath').text(inputs.inputDwellPath);
	$('#lblMatchVelocity').text(inputs.matchVelocity == 1 ? 'Y' : 'N');
	$('#lblDwellCatchScale').text(inputs.dwellCatchScale);
	$('#lblDwellTossScale').text(inputs.dwellTossScale);
	$('#lblEmptyTossScale').text(inputs.emptyTossScale);
	$('#lblEmptyCatchScale').text(inputs.emptyCatchScale);
	$('#lblArmAngle').text(inputs.armAngle + ' rad');

	var lblSurfaces = "";
	for (var i = 0; i < inputs.surfaces.length; i++) {
		var surface = inputs.surfaces[i];
		lblSurfaces += 'surface ' + i + ': position <' + surface.position.x + ',' + surface.position.y + ',' + surface.position.z + '>' + ' normal <' + surface.normal.x + ',' + surface.normal.y + ',' + surface.normal.z + '> half-width ' + surface.scale + 'm';
		if (i < inputs.surfaces.length-1) {
			lblSurfaces += ", ";
		}
	}
	$('#lblSurfaces').text(lblSurfaces);

}

function runExample(exampleName) {
	$.getJSON("examples.json", function(data) {
		for (var i = 0; i < data.examples.length; i++) {
			if (data.examples[i].name == exampleName) {
				bindInputs(applyInputDefaults(data.examples[i]));
				go();
			}
		}
	});
}

function buildExamples() {
	$.getJSON("examples.json", function(data) {
		for (var i = 0; i < data.examples.length; i++) {
			$('#examplesList').append('<li><a href="#" onclick="runExample(\'' + data.examples[i].name + '\');">' + data.examples[i].name + '</a></li>');
		}
	});
}

function generateGIF() {

	$('#gifProgress').show();
	$('#gifLink').empty();

	animator.paused = true;
	var numFrames = Math.round((siteswap.states.length*siteswap.beatDuration)*35);
	var currentFrame = 0;

	var canvas = document.createElement( 'canvas' );
	canvas.width = animator.renderer.domElement.width;
	canvas.height = animator.renderer.domElement.height;

	var context = canvas.getContext( '2d' );

	var buffer = new Uint8Array( canvas.width * canvas.height * numFrames * 5 );
	var gif = new GifWriter( buffer, canvas.width, canvas.height, { loop: 0 } );

	var pixels = new Uint8Array( canvas.width * canvas.height );

	var addFrame = function () {

		animator.render((currentFrame/numFrames)*(siteswap.states.length*siteswap.beatDuration)*1000);

		context.drawImage( animator.renderer.domElement, 0, 0 );

		var data = context.getImageData( 0, 0, canvas.width, canvas.height ).data;

		var palette = [];

		for ( var j = 0, k = 0, jl = data.length; j < jl; j += 4, k ++ ) {

			var r = Math.floor( data[ j + 0 ] * 0.1 ) * 10;
			var g = Math.floor( data[ j + 1 ] * 0.1 ) * 10;
			var b = Math.floor( data[ j + 2 ] * 0.1 ) * 10;
			var color = r << 16 | g << 8 | b << 0;

			var index = palette.indexOf( color );

			if ( index === -1 ) {

				pixels[ k ] = palette.length;
				palette.push( color );

			} else {

				pixels[ k ] = index;

			}

		}

		// force palette to be power of 2

		var powof2 = 1;
		while ( powof2 < palette.length ) powof2 <<= 1;
		palette.length = powof2;

		gif.addFrame( 0, 0, canvas.width, canvas.height, pixels, { palette: new Uint32Array( palette ), delay: 5 } );

		$('#gifProgress').val(currentFrame/numFrames);

		currentFrame++;
		if (currentFrame == numFrames) {
			finish();
		} else {
			setTimeout(addFrame,0);
		}


	}

	var finish = function () {

		// return buffer.slice( 0, gif.end() );

		var string = '';

		for ( var i = 0, l = gif.end(); i < l; i ++ ) {

			string += String.fromCharCode( buffer[ i ] )

		}

		$('#gifLink').append("<a href='" + 'data:image/gif;base64,' + btoa( string ) + "' target='_blank'>Download GIF</a>");

		animator.paused = false;
		animator.animate();
		$('#gifProgress').hide();

	}

	addFrame();

}

function findSiteswaps() {
	window.onresize();
	$('#siteswapsList').empty();
	$('#graphContainer').empty();	
	$('#activeSiteswapContainer').hide();

	var excludeInput = $('#exclude').val();

	var config = {
		minPeriod: 1, 
		maxPeriod: parseInt($("#explorerMaxPeriod").val()),
		numProps: parseInt($("#explorerNumProps").val()),
		maxSiteswaps: parseInt($("#explorerMaxSiteswaps").val()),
		includeExcited: $('#explorerIncludeExcited')[0].checked,
		includeMultiplex: $('#explorerIncludeMultiplex')[0].checked,
		async: true,
		sync: $('#explorerSync')[0].checked,
		callbacks: {
			siteswapFound: function (siteswap, siteswapIx, excited) {
				$('#siteswapsList').append('<li><a class="' + (excited ? 'excited' : 'ground') + '" href="#" onclick="runSiteswap(\''+siteswap+'\')">'+siteswap+'</a></li>');
			}
		}
	};
	
	SiteswapGraph.siteswapGraph(config);
}

function runSiteswap(s) {
	$('#siteswap').val(s);
	updateAdvancedInputsFromBasic();
	go();
}

function updateDrawHandsForProp() {
	if ($('#prop').val() != 'ball') {
		$('#drawHands')[0].checked = false;
	}
}

function showHideCameraCustomPosition() {
	if ($('#cameraMode').val() == "custom") {
		$('#cameraCustomPositionContainer').show();
	} else {
		$('#cameraCustomPositionContainer').hide();
	}
}

function getInputsFromQueryString() {
	var inputs = {};
	var siteswap = util.getURLQueryStringParameterByName("siteswap");
	var props = JSON.parse(util.getURLQueryStringParameterByName("props"));
	var beatDuration = util.getURLQueryStringParameterByName("beatDuration");
	var dwellPath = util.getURLQueryStringParameterByName("dwellPath");
	
	if(siteswap !== null) {
		inputs.siteswap = siteswap;
	}
	if (props !== null) {
		inputs.props = props;
	}
	if (beatDuration !== null) {
		inputs.beatDuration = beatDuration;
	}
	if (dwellPath !== null) {
		inputs.dwellPath = dwellPath;
	}

	return inputs;
}

function saveCurrentSiteswap() {
	var savedSiteswaps = getSavedSiteswaps();
	savedSiteswaps.push({name: $('#savedName').val(), version: 16, inputs: $('#inputsAdvanced').val()});
	window.localStorage.setItem("savedSiteswaps",JSON.stringify(savedSiteswaps));
	refreshSavedSiteswapsList();
	$('#saveSuccess').show(1000);
	setTimeout(function() { $('#saveSuccess').hide(1000); }, 2000);
}

function getSavedSiteswaps() {	
	if (window.localStorage.getItem("savedSiteswaps") === null) {
		window.localStorage.setItem("savedSiteswaps", "[]");
	} 
	return JSON.parse(window.localStorage.getItem("savedSiteswaps"));	
}

function refreshSavedSiteswapsList() {
	var savedSiteswaps = getSavedSiteswaps();
	var savedList = $('#savedList');
	savedList.empty();
	for(var i = 0; i < savedSiteswaps.length; i++) {
		savedList.append('<li><a href="#" class="deleteLink" onclick="deleteSavedSiteswap(' + i + ');"><span class="glyphicon glyphicon-remove"></span></a><a href="#" onclick="runSavedSiteswap(' + i + ');">'+ savedSiteswaps[i].name +'</a></li>');
	}	
}

function runSavedSiteswap(savedSiteswapIndex) {
	var savedSiteswaps = getSavedSiteswaps();
	$('#inputsAdvanced').val(savedSiteswaps[savedSiteswapIndex].inputs);
	go();
}

function deleteSavedSiteswap(savedSiteswapIndex) {
	var savedSiteswaps = getSavedSiteswaps();
	savedSiteswaps.splice(savedSiteswapIndex);
	window.localStorage.setItem("savedSiteswaps",JSON.stringify(savedSiteswaps));
	refreshSavedSiteswapsList();
}

window.go = go;
window.displayMenu = displayMenu;
window.updateAdvancedInputsFromBasic = updateAdvancedInputsFromBasic;
window.zoomIn = zoomIn;
window.zoomOut = zoomOut;
window.updateAnimationSpeed = updateAnimationSpeed;
window.updateCameraMode = updateCameraMode;
window.updateDisplayPropPaths = updateDisplayPropPaths;
window.updateAdvancedLabels = updateAdvancedLabels;
window.runExample = runExample;
window.generateGIF = generateGIF;
window.findSiteswaps = findSiteswaps;
window.runSiteswap = runSiteswap;
window.updateDrawHandsForProp = updateDrawHandsForProp;
window.showHideCameraCustomPosition = showHideCameraCustomPosition;
window.saveCurrentSiteswap = saveCurrentSiteswap;
window.runSavedSiteswap = runSavedSiteswap;
window.deleteSavedSiteswap = deleteSavedSiteswap;
},{"./Siteswap.js":3,"./SiteswapAnimator.js":4,"./SiteswapGraph.js":5,"./util.js":7}],7:[function(require,module,exports){
module.exports.getURLQueryStringParameterByName = function(name) {
    name = name.replace(/[\[]/, "\\[").replace(/[\]]/, "\\]");
    var regex = new RegExp("[\\?&]" + name + "=([^&#]*)"),
        results = regex.exec(location.search);
    return results === null ? null : decodeURIComponent(results[1].replace(/\+/g, " "));
}

/* used for deep cloning of various arrays/objects */
module.exports.cloneObject = function(obj) {
  var newObj = (obj instanceof Array) ? [] : {};
  for (var i in obj) {
    if (i == 'clone') continue;
    if (obj[i] && typeof obj[i] == "object") {
      newObj[i] = module.exports.cloneObject(obj[i]);
    } else newObj[i] = obj[i]
  } return newObj;
}

/* used to check if two states are the same */
module.exports.arraysEqual = function(a,b) {
	if (a instanceof Array && b instanceof Array) {
		if (a.length != b.length) {
			return false;
		} else {
			for (var i = 0; i < a.length; i++) {
				/* if this is a multi-dimensional array, check equality at the next level */
				if (a[i] instanceof Array || b[i] instanceof Array) {
					var tmp = module.exports.arraysEqual(a[i],b[i]);
					if (!tmp) {
						return false;
					}
				} else if (a[i] != b[i]) {
					return false;
				}
			}
		}
	} else {
		return false;
	}
	return true;
}

if (!Array.prototype.last){
    Array.prototype.last = function(){
        return this[this.length - 1];
    };
};

module.exports.cross = function(a,b) {
	return {
		x: a.y*b.z - a.z*b.y,
		y: a.z*b.x - a.x*b.z,
		z: a.x*b.y - a.y*b.x
	}
}

module.exports.dot = function(a,b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

module.exports.magnitude = function(a) {
	return Math.sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

module.exports.normalize = function(a) {
	var mag = module.exports.magnitude(a);
	a.x = a.x/mag;
	a.y = a.y/mag;
	a.z = a.z/mag;
	return a;
}

module.exports.multiply = function(a,b) {
	a.x *= b;
	a.y *= b;
	a.z *= b;
	return a;
}

module.exports.add = function(a,b) {
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}

module.exports.sub = function(a,b) {
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}

module.exports.negate = function(a) {
	module.exports.multiply(a,-1);
	return a;
}
},{}]},{},[6]);
