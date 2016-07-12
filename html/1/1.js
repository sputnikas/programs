function showMessage(){
	alert ("Вы щелкнули по div-у");
}

function areaRectangle(obj){
	var a = obj.t1.value;
	var b = obj.t2.value;
	var s = a*b;
	obj.res.value = s;
}

function message(m){
	alert (m);
}

function showDesc(obj, n){
	obj.desc.value=n;
}

function delet(obj){
	obj.desc.value='';
}

function areaOfTriangle(obj){
	var a=1*obj.st1.value;
	var b=1*obj.st2.value;
	var c=1*obj.st3.value;
	var p=(a+b+c)/2;  
	var s=Math.sqrt(p*(p-a)*(p-b)*(p-c));
	obj.res.value=s;
}

function resize(obj){
	var w=obj.width;
	if (w<1000) {
		obj.width*=2;
		setTimeout("bigPict()", 500)
	}
}

function repair(obj){
	var w=obj.width;
	if (w<1000) {
		obj.width/=2;
		setTimeout("bigPict()", 500)
	}
}