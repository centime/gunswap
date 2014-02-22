/* --------- */
/* CONSTANTS */
/* --------- */
var LEFT = 0, RIGHT = 1;

/* ------- */
/* HELPERS */
/* ------- */

Object.prototype.clone = function() {
  var newObj = (this instanceof Array) ? [] : {};
  for (i in this) {
    if (i == 'clone') continue;
    if (this[i] && typeof this[i] == "object") {
      newObj[i] = this[i].clone();
    } else newObj[i] = this[i]
  } return newObj;
};

/* 
	helper for figuring out the number of props 
	TODO : this should interpret a-g as valid toss heights
*/
function sumIntegers(str) {
	if (str.length == 1 && parseInt(str)) {
		return parseInt(str);
	} else if (str.length > 1) {
		return str.split('').reduce(function(prev,cur) { return (parseInt(prev) ? parseInt(prev) : 0) + (parseInt(cur) ? parseInt(cur) : 0) });
	} else {
		return 0;
	}
}

/* helper for state array equality */
function arraysEqual(a,b) {
	if (a instanceof Array && b instanceof Array) {
		if (a.length != b.length) {
			return false;
		} else {
			for (var i = 0; i < a.length; i++) {
				/* if this is a multi-dimensional array, check equality at the next level */
				if (a[i] instanceof Array || b[i] instanceof Array) {
					var tmp = arraysEqual(a[i],b[i]);
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

/* returns a string describing a given state */
function printState(state) {

	var str = '';

	for(var juggler = 0; juggler < state.length; juggler++) {
		str += ('J' + juggler);
		for (var hand = 0; hand < state[juggler].length; hand++) {
			str += (hand == LEFT ? ' L ' : ' R ');
			for (var beat = 0; beat < state[juggler][hand].length; beat++) {
				str += (state[juggler][hand][beat] == undefined ? 'X' : state[juggler][hand][beat][0]);
			}
		}
		str += ' ';
	}

	return str;
}

/* ---------------*/
/* SITESWAP CLASS */
/* -------------- */

function Siteswap(siteswap) {	

	this.siteswap = siteswap;

	/* first check if the siteswap is a passing siteswap, that will inform the pattern for a valid toss */
	this.numJugglers = 1;
	var isPassingPattern = /<[^ ]+>/.test(siteswap);

	var numJugglerMismatch = false;

	if (isPassingPattern) {
		var passingBeatArray = siteswap.match(/<[^ <>]+>/g);
		this.numJugglers = passingBeatArray[0].split("|").length;

		/* 
			check to make sure each beat in the passing pattern has the same number of jugglers 
			if a passing pattern only has 1 juggler than it's automatically a mismatch
		*/
		numJugglerMismatch = (this.numJugglers == 1 ? true : false);
		
		var numJugglersTmp = this.numJugglers;
		passingBeatArray.map(function(a) { 
			if (a.split("|").length != numJugglersTmp) 
				{ numJugglerMismatch = true; } 
		});
	}

	/* the number of jugglers determines a valid pass pattern */
	var passPattern = "";
	if (this.numJugglers == 2) {
		passPattern = "p";
	} else if (this.numJugglers > 2) {
		passPattern = "p[1-" + this.numJugglers + "]";
	}

	/* construct the various regex patterns. see blog post for details about this */
	var validToss = "(\\d|[a-g])x?(" + passPattern + ")?";
	var validMultiplex = "\\[(" + validToss + ")+\\]";
	var validSync = "\\((" + validToss + "|" + validMultiplex + "),(" + validToss + "|" + validMultiplex + ")\\)";
	var validBeat = "(" + validToss + "|" + validMultiplex + "|" + validSync + ")";
	var validPass = "<" + validBeat + "(\\|" + validBeat + ")+>";
	var validSiteswap = "^(" + validPass + ")+|(" + validBeat + ")+$";

	var validTossRe = new RegExp(validToss,"g");
	var validMultiplexRe = new RegExp(validMultiplex,"g");
	var validSyncRe = new RegExp(validSync,"g");
	var validBeatRe = new RegExp(validBeat,"g");
	var validPassRe = new RegExp(validPass,"g");
	var validSiteswapRe = new RegExp(validSiteswap,"g");	

	if (siteswap.match(validSiteswapRe) && !numJugglerMismatch) {
		
		/* get the array of each beats' tosses */
		var beatArr = isPassingPattern ? siteswap.match(validPassRe) : siteswap.match(validBeatRe);
		
		/* figure out how many props */
		var tmp = 0;
		beatArr.map(function(beat) {
			if (beat.match(validPassRe)) {
				var patterns = beat.split('|');
				for (var i = 0; i < patterns.length; i++) {
					if (i == 0) {
						patterns[i] = patterns[i].substr(1);
					} 
					if (i == patterns.length-1) {
						patterns[i] = patterns[i].substr(0,patterns[i].length-1);
					}
					if (patterns[i].match(validSyncRe)) {
						tmp += sumIntegers(patterns[i])/2;
					} else {
						tmp += sumIntegers(patterns[i]);
					}
				}
			} else {
				if (beat.match(validSyncRe)) {
					tmp += sumIntegers(beat)/2;
				} else {
					tmp += sumIntegers(beat);
				}				
			}
		});

		if((tmp/beatArr.length % 1) == 0) {
			this.numProps = tmp/beatArr.length;
		}

		/* if there is a valid number of props, start creating the state array */
		if(this.numProps) {

			/* initialize array of tosses */
			this.tossArr = [];

			/* for each beat get all the tosses */
			for (var i = 0; i < beatArr.length; i++) {
				var tosses = [];
				getTosses(tosses,beatArr[i], 0 /* assume juggler 0 */);
				this.tossArr.push(tosses);
			}

			/* figure out the max throw height which will inform the size of the state array */
			this.maxTossHeight = 0;

			for (var i = 0; i < this.tossArr.length; i++) {
				for (var j = 0; j < this.tossArr[i].length; j++) {
					if(this.tossArr[i][j].numBeats > this.maxTossHeight) {
						this.maxTossHeight = this.tossArr[i][j].numBeats;
					}
				}
			}			
			
			/* create a queue of props */
			var props = [];

			for (var i = 0; i < this.numProps; i++) {
				props.push(i);
			}

			/* initialize the state array */
			this.states = [];

			/* initialize current state */
			var curState = [];
			for (var j = 0; j < this.numJugglers; j++) {
				curState.push([[],[]]);
				for (var k = 0; k < this.maxTossHeight; k++) {
					curState[j][LEFT].push(undefined);
					curState[j][RIGHT].push(undefined);
				}
			}

			var patternComplete = false;
			var initComplete = false;
			var beat = 0;
			var hand = LEFT;

			/* keep going until pattern complete */
			while (!patternComplete) {

				/* queue of props to throw this beat */
				var propsLanding = [];

				/* update the current state for each juggler */
				for (var j = 0; j < this.numJugglers; j++) {
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
				for (var j = 0; j < this.tossArr[beat % this.tossArr.length].length; j++) {
					
					var toss = this.tossArr[beat % this.tossArr.length][j];
					var tossHand = (toss.hand == undefined ? hand : toss.hand);
					var catchHand = (toss.crossing ? 1 - tossHand : tossHand);

					var prop = undefined;

					/* iterate through the props landing and look for one landing in the hand that this toss is occurring */
					for (var k = 0; k < propsLanding.length; k++) {
						if(propsLanding[k].juggler == toss.juggler && propsLanding[k].hand == tossHand) {
							prop = propsLanding[k].propId;
						}
					}

					/* if no props landing to be thrown */
					if (prop == undefined) {
						prop = props.shift();			
					} 

					/* if prop is still undefined (ie. there are none left) then we've got an invalid siteswap */
					if (prop == undefined) {
						throw 'Invalid';
					}

					if(curState[toss.targetJuggler][catchHand][toss.numBeats-1] == undefined) {
						curState[toss.targetJuggler][catchHand][toss.numBeats-1] = [prop];
					} else {
						curState[toss.targetJuggler][catchHand][toss.numBeats-1].push(prop);
					}
				}
								

				/* if we're at the beginning of the toss array and we've returned to the original state, the pattern is complete */
				if (initComplete && beat % this.tossArr.length == 0 && arraysEqual(this.states[0],curState)) {					
					patternComplete = true;				
				} else {
					/* add the current state to the state array */
					this.states.push(curState.clone());
				}					

				/* if all props have been introduced to pattern and we're at the end of the pattern, init is complete and steady-state pattern truly begins with the next beat */
				if (props.length == 0 && (beat+1) % this.tossArr.length == 0 && !initComplete) {
					initComplete = true;
					this.states = []; /* reset the state */
				}			

				beat++;
				hand = 1 - hand; //alternate hands

				/* fail safe in case the pattern is too long */
				if (beat > 100) {
					patternComplete = true;
				}

			}

		}

	} else {
		throw 'Invalid siteswap format';
	}	

	/* get all the tosses for a given beat's siteswap */
	function getTosses(tosses, siteswap, juggler, sync, hand) {
		if (siteswap.match(validPassRe)) {
			var patterns = siteswap.match(validBeatRe);
			patterns.map(function(s,ix) {
				getTosses(tosses, s, ix);				
			});
		} else if (siteswap.match(validSyncRe)) {
			var patterns = siteswap.split(",");
			getTosses(tosses,patterns[0].substr(1),juggler,true,LEFT);
			getTosses(tosses,patterns[1].substr(0,patterns.length),juggler,true,RIGHT);
		} else if (siteswap.match(validMultiplex)) {
			var patterns = siteswap.match(validTossRe);
			patterns.map(function(s) {
				getTosses(tosses,s,juggler);
			});
		} else {
			var numBeats = parseInt(siteswap[0]);
			var targetJuggler = juggler;

			var pIx = siteswap.indexOf("p");
			if (pIx > 0) {
				if (this.numJugglers > 2) {
					targetJuggler = siteswap[pIx+1];
				} else {
					targetJuggler = 1 - juggler;
				}			
			}

			var crossing = (siteswap.length > 1 && (siteswap[1] == "x" && numBeats % 2 == 0) || (siteswap[1] != "x" && numBeats % 2 == 1)) || numBeats % 2 == 1 ? true : false;

			if (sync) {
				numBeats = numBeats/2;
			}

			tosses.push(
				{
					juggler: juggler,
					targetJuggler: targetJuggler,
					hand: hand, /* only !undefined for sync throws */
					crossing: crossing,
					numBeats: numBeats,
					siteswap: siteswap
				}
			);
		}
	}

	this.debugStates = function() {
		for (var i = 0; i < this.states.length; i++) {
			$('#states').append(printState(this.states[i]) + '</br>');
		}
	}

}