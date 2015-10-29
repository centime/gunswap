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