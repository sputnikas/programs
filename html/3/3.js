var canvas;
var context;
var timer;

window.onload = function() {
	canvas = document.getElementById("canvas");
	canvas.height = 454;
	canvas.width = 600;
	context = canvas.getContext("2d");
	timer = -1;
}

function Particle(m, x, y, vx, vy, R) {
	this.m = m || 0;
	this.x = x || 0;
	this.y = y || 0;
	this.vx = vx || 0;
	this.vy = vy || 0;
	this.R = R || 0;
}

var p = [];
var dt;

function beginDrawing(obj) {
	p = [];
	var m = 1.*obj.m1.value;
	var x = 1.*obj.x1.value;
	var y = 1.*obj.y1.value;
	var vx = 1.*obj.vx1.value;
	var vy = 1.*obj.vy1.value;
	var R = 1.*obj.R1.value;
	p.push(new Particle(m, x, y, vx, vy, R));
	m = 1.*obj.m2.value;
	x = 1.*obj.x2.value;
	y = 1.*obj.y2.value;
	vx = 1.*obj.vx2.value;
	vy = 1.*obj.vy2.value;
	R = 1.*obj.R2.value;
	p.push(new Particle(m, x, y, vx, vy, R));
	dt = 1.*obj.dt.value/100;
	if (timer != -1) clearTimeOut(timer);
	timer = setTimeout("drawFrame()", dt*20);
	obj.button3.style.display = "none";
	obj.button2.style.display = "inline";
}

function stopDrawing(obj) {
	if (timer != -1) clearTimeout(timer);
	timer = -1;
	obj.button2.style.display = "none";
	obj.button3.style.display = "inline";
}

function continueDrawing(obj) {
	timer = setTimeout("drawFrame()", dt*20);
	obj.button3.style.display = "none";
	obj.button2.style.display = "inline";
}

function drawFrame() {
	context.clearRect(0, 0, canvas.width, canvas.height);
	context.beginPath();
	
	for (var i=0; i<p.length; i++) {
		context.arc(p[i].x, p[i].y, p[i].R, 0, Math.PI*2);
		context.lineStyle = "#109bfc";
		context.lineWidth = 1;
		context.fillStyle = "red";
		context.fill();
	}
	
	var rx, ry, vx, vy, v, r;
	var dx, dy, t, ct, R, dpx, dpy, nx, ny;
	var energy = 0.;
	
	for (var i=0; i<p.length; i++) {
		for (var j=0; j<p.length; j++) {
			if (j != i) {
				rx = p[i].x - p[j].x;
				ry = p[i].y - p[j].y;
				vx = p[i].vx - p[j].vx;
				vy = p[i].vy - p[j].vy;
				dx = rx + vx*dt;
				dy = ry + vy*dt;
				R = p[i].R + p[j].R;
				r = Math.sqrt(Math.pow(rx, 2) + Math.pow(ry, 2));
				v = Math.sqrt(Math.pow(vx, 2) + Math.pow(vy, 2));
				ct = (rx*vx + ry*vy)/r/v;
				if ((ct < 0) && (r/R > 1) && (R/r > Math.sqrt(1 - ct*ct))) {
					t = - (ct + Math.sqrt(R*R/r/r - 1 + ct*ct))*r/v;
					if (t < dt) {
						nx = (rx + vx*t)/R;
						ny = (ry + vy*t)/R;
						dpx = 2.*p[i].m*p[j].m/(p[i].m + p[j].m)*(vx*nx + vy*ny)*nx;
						dpy = 2.*p[i].m*p[j].m/(p[i].m + p[j].m)*(vx*nx + vy*ny)*ny;
						p[i].vx -= dpx/p[i].m;
						p[j].vx += dpx/p[j].m;
						p[i].vy -= dpy/p[i].m;
						p[j].vy += dpy/p[j].m;
					}
				}
			}
		}
		if ((p[i].x + p[i].vx*dt > canvas.width - p[i].R) || (p[i].x + p[i].vx*dt < p[i].R)) {
			if (p[i].x + p[i].vx*dt > canvas.width - p[i].R) { 
				p[i].x = canvas.width - p[i].R - (dt - (canvas.width - p[i].R - p[i].x)/p[i].vx)*p[i].vx; 
			} else {
				p[i].x = p[i].R - (dt + (p[i].x - p[i].R)/p[i].vx)*p[i].vx; 
			}
			p[i].vx *= -1;
		} else {
			p[i].x = p[i].x + p[i].vx*dt;
		}
		if ((p[i].y + p[i].vy*dt > canvas.height - p[i].R) || (p[i].y + p[i].vy*dt < p[i].R)) {
			if (p[i].y + p[i].vy*dt > canvas.height - p[i].R) { 
				p[i].y = canvas.height - p[i].R - (dt - (canvas.height - p[i].R - p[i].y)/p[i].vy)*p[i].vy; 
			} else {
				p[i].y = p[i].R - (dt + (p[i].y - p[i].R)/p[i].vy)*p[i].vy; 
			}
			p[i].vy *= -1;
		} else {
			p[i].y = p[i].y + p[i].vy*dt;
		}
		energy += p[i].m*(p[i].vx*p[i].vx + p[i].vy*p[i].vy)/2;
	}
	
	data1.result.value = energy;
	
	timer = setTimeout("drawFrame()", dt*20);
}