window.onload = function () {

	displayMenu('Examples');

	updateAdvancedInputsFromBasic();
	updateAdvancedLabels();

	$('#menuContainer').height($(window).height());

	window.animator = new SiteswapAnimator.SiteswapAnimator('animatorCanvasContainer');

	buildExamples();

}

window.onresize = function () {
	//animator.resize($('#animatorContainer').width()-25, $(window).height()-40);
	$('#menuContainer').height($(window).height());
}

function displayMenu(menu) {
	$('#menuBasic').hide();
	$('#menuAdvanced').hide();
	$('#menuExamples').hide();
	$('#menuHelp').hide();
	$('#menuAbout').hide();
	$('#menu' + menu).show();

	$('#navBasic').removeClass('activeNav');
	$('#navAdvanced').removeClass('activeNav');
	$('#navExamples').removeClass('activeNav');
	$('#navHelp').removeClass('activeNav');
	$('#navAbout').removeClass('activeNav');
	$('#nav' + menu).addClass('activeNav');	

	if (menu == 'Basic') {
		updateAdvancedInputsFromBasic();
	}

}

function updateAdvancedInputsFromBasic() {
	bindInputs({
		siteswap: $('#siteswap').val(),
		props: [{type: $('#prop').val(), color: 'red', radius: .05, C: .97}],
		beatDuration: $('#beatDuration').val(),
		dwellRatio: .5,
		dwellPath: $('#dwellPath').val(),
		matchVelocity: 1,
		dwellCatchScale: .06,
		dwellTossScale: .06,
		emptyTossScale: .025,
		emptyCatchScale: .025,
		armAngle: .1,
		surfaces: [
			{
				position: {x:0, y:0, z:0},
				normal: {x:0, y:1, z:0},
				scale: 5,
			}
		]
	});
}

function applyInputDefaults(inputs) {
	inputs.siteswap = inputs.siteswap === undefined ? "3" : inputs.siteswap;
	inputs.props = inputs.props === undefined ? [{type: "ball", color: "red", radius: ".05", C: .97}] : inputs.props;
	inputs.beatDuration = inputs.beatDuration === undefined ? .25 : inputs.beatDuration;
	inputs.dwellRatio = inputs.dwellRatio === undefined ? .5 : inputs.dwellRatio;
	inputs.dwellPath = inputs.dwellPath === undefined ? "(30)(10)" : inputs.dwellPath;
	inputs.matchVelocity = inputs.matchVelocity === undefined ? 1 : inputs.matchVelocity;
	inputs.dwellCatchScale = inputs.dwellCatchScale === undefined ? .06 : inputs.dwellCatchScale;
	inputs.dwellTossScale = inputs.dwellTossScale === undefined ? .06 : inputs.dwellTossScale;
	inputs.emptyTossScale = inputs.emptyTossScale === undefined ? .025 : inputs.emptyTossScale;
	inputs.emptyCatchScale = inputs.emptyCatchScale === undefined ? .025 : inputs.emptyCatchScale;
	inputs.armAngle = inputs.armAngle === undefined ? .1 : inputs.armAngle;
	inputs.surfaces = inputs.surfaces === undefined ? [{position: {x:0, y:0, z:0}, normal: {x:0, y:1, z:0}, scale: 5}] : inputs.surfaces;
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
	for (var i = 0; i < inputs.surfaces.length; i++) {
		inputsText += inputs.surfaces[i].position.x + " " + inputs.surfaces[i].position.y + " " + inputs.surfaces[i].position.z + " " + inputs.surfaces[i].normal.x + " " + inputs.surfaces[i].normal.y + " " + inputs.surfaces[i].normal.z + " " + inputs.surfaces[i].scale + "\n";
	}
	$('#inputsAdvanced').val(inputsText);
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

	var surfaces= [];
	for (var i = 5; i < lines.length; i++) {
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
		surfaces: surfaces
	};
}

function go() {

	var inputs = parseInputs($('#inputsAdvanced').val());

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
			surfaces: 			inputs.surfaces
		});

	animator.go(siteswap, {});
}

function zoomIn() { animator.zoomIn(); }

function zoomOut() { animator.zoomOut(); }

function updateAnimationSpeed() {
	animationSpeed = parseFloat($('#animationSpeed').val());
	animator.updateAnimationSpeed(animationSpeed);
}

function updateCameraMode() {
	cameraMode = $('#cameraMode').val();
	animator.updateCameraMode(cameraMode);
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
			$('#examples').append('<li><a href="#" onclick="runExample(\'' + data.examples[i].name + '\');">' + data.examples[i].name + '</a></li>');
		}
	});
}

function filterExamples() {
	var filterValue = $('#filterExamples').val();
	if (filterValue == '') {
		$('#examples li').show();	
	} else {
		$('#examples li').hide();
		$('#examples li:contains("' + filterValue + '")').show();
	}
}