
(function(){
"use strict";
var $=function(id){return document.getElementById(id)};
var timer=null,lastBundle=null;
var colors={onoff:"#dc2626",p:"#07855b",pi:"#7c3aed",pid:"#2563eb",base:"#64748b",sp:"#6b7280",output:"#f59e0b"};
var historyColors=["#0f766e","#9333ea","#c2410c","#0369a1","#a21caf","#4d7c0f"];
var histories={p:[],pi:[],pid:[]};
var defaults={
 model:"fopdt",initial:25,bias:25,tau:300,gain:16,delay:30,cinitial:25,cdelay:30,equation:"(-(y - 25) + 16*u) / 300",
 sp:1000,duration:5400,dt:1,umin:0,umax:100,hys:5,pk:0.3,pikc:0.3,piti:200,pidkc:0.3,pidti:200,pidtd:25,aw:"on"
};
var fixed={p:{kp:0.3},pi:{kc:0.3,ti:200},pid:{kc:0.3,ti:200,td:25}};

function num(id,fb){var v=Number($(id).value);return Number.isFinite(v)?v:fb}
function clamp(v,a,b){return Math.min(b,Math.max(a,v))}
function fmt(v,d){if(v==null||!Number.isFinite(v))return "—";return v.toLocaleString("ko-KR",{minimumFractionDigits:d,maximumFractionDigits:d})}
function val(id,v){$(id).value=v}

function getConfig(){
 var c={
  model:$("model").value,sp:num("sp",1000),duration:num("duration",5400),dt:num("dt",1),
  umin:num("umin",0),umax:num("umax",100),hys:Math.max(0,num("hys",5)),aw:$("aw").value==="on"
 };
 if(!(c.dt>0)||!(c.duration>0)||c.duration/c.dt>120000)throw Error("dt 또는 총 시간을 확인하세요. 최대 120,000 step입니다.");
 if(!(c.umax>c.umin))throw Error("출력 최대값은 최소값보다 커야 합니다.");
 if(c.model==="fopdt"){
  c.initial=num("initial",25);c.bias=num("bias",25);c.tau=num("tau",300);c.gain=num("gain",16);c.delay=Math.max(0,num("delay",30));
  if(!(c.tau>0))throw Error("시정수 τ는 0보다 커야 합니다.");
 }else{
  c.initial=num("cinitial",25);c.delay=Math.max(0,num("cdelay",30));c.eq=$("equation").value.trim();
  if(!c.eq)throw Error("Custom ODE 식을 입력하세요.");
  if(!/^[0-9A-Za-z_+\-*/()., ^]+$/.test(c.eq))throw Error("Custom ODE에 허용되지 않은 문자가 있습니다.");
  c.fn=new Function("y","u","t","sp","Math","return ("+c.eq+")");
 }
 return c;
}

function currentParams(kind){
 if(kind==="p")return{kp:Math.max(0,num("pkn",0.3))};
 if(kind==="pi")return{kc:Math.max(0,num("pikcn",0.3)),ti:Math.max(0.0001,num("pitin",200))};
 return{kc:Math.max(0,num("pidkcn",0.3)),ti:Math.max(0.0001,num("pidtin",200)),td:Math.max(0,num("pidtdn",25))};
}

function derivative(c,y,u,t){
 if(c.model==="fopdt")return (-(y-c.bias)+c.gain*u)/c.tau;
 var z=c.fn(y,u,t,c.sp,Math);
 if(!Number.isFinite(z))throw Error("Custom ODE 계산 결과가 유한수가 아닙니다.");
 return z;
}

function simulate(c,kind,p){
 var steps=Math.floor(c.duration/c.dt)+1;
 var t=new Float64Array(steps),y=new Float64Array(steps),u=new Float64Array(steps),e=new Float64Array(steps);
 y[0]=c.initial;
 var delaySteps=Math.max(0,Math.round(c.delay/c.dt));
 var integ=0,prevErr=c.sp-y[0],prevU=c.umax,satTime=0;
 for(var i=0;i<steps;i++){
  t[i]=i*c.dt;
  var err=c.sp-y[i],ctrl=0;
  e[i]=err;
  if(kind==="onoff"){
   if(y[i]<c.sp-c.hys)ctrl=c.umax;
   else if(y[i]>c.sp+c.hys)ctrl=c.umin;
   else ctrl=prevU;
  }else if(kind==="p"){
   ctrl=clamp(p.kp*err,c.umin,c.umax);
  }else{
   var de=i?(err-prevErr)/c.dt:0;
   var cand=integ+err*c.dt;
   var td=kind==="pi"?0:p.td;
   var raw=p.kc*(err+cand/p.ti+td*de);
   ctrl=clamp(raw,c.umin,c.umax);
   var saturated=Math.abs(ctrl-raw)>1e-10;
   if(!(c.aw&&saturated))integ=cand;
   if(c.aw&&saturated){
    raw=p.kc*(err+integ/p.ti+td*de);
    ctrl=clamp(raw,c.umin,c.umax);
   }
  }
  u[i]=ctrl;prevU=ctrl;
  if(ctrl<=c.umin+1e-9||ctrl>=c.umax-1e-9)satTime+=c.dt;
  if(i<steps-1){
   var j=i-delaySteps,delayedU=j>=0?u[j]:0;
   var next=y[i]+c.dt*derivative(c,y[i],delayedU,t[i]);
   if(!Number.isFinite(next)||Math.abs(next)>1e12)throw Error("응답이 발산했습니다. 공정식 또는 제어기 값을 확인하세요.");
   y[i+1]=next;
  }
  prevErr=err;
 }
 return{t:t,y:y,u:u,e:e,sat:satTime,kind:kind};
}

function perf(c,r){
 var max=-Infinity,min=Infinity,iae=0;
 for(var i=0;i<r.y.length;i++){max=Math.max(max,r.y[i]);min=Math.min(min,r.y[i]);iae+=Math.abs(r.e[i])*c.dt}
 var fin=r.y[r.y.length-1],err=c.sp-fin,span=Math.abs(c.sp-c.initial),ov=span?Math.max(0,(max-c.sp)/span*100):0;
 var band=Math.max(Math.abs(c.sp)*0.02,1e-9),last=-1;
 for(i=0;i<r.y.length;i++)if(Math.abs(r.y[i]-c.sp)>band)last=i;
 var settling=last<r.y.length-1?r.t[Math.min(last+1,r.t.length-1)]:null;
 return{fin:fin,err:err,max:max,min:min,iae:iae,ov:ov,settling:settling,sat:r.sat};
}

function lastWindowWidth(c,r,seconds){
 var start=Math.max(0,c.duration-seconds),mn=Infinity,mx=-Infinity;
 for(var i=0;i<r.t.length;i++)if(r.t[i]>=start){mn=Math.min(mn,r.y[i]);mx=Math.max(mx,r.y[i])}
 return{width:mx-mn,min:mn,max:mx};
}

function canvasSetup(cv,heightRef){
 var w=Math.max(320,cv.clientWidth||800),ratio=Math.max(1,window.devicePixelRatio||1);
 var h=w*heightRef/1200;
 cv.width=Math.round(w*ratio);cv.height=Math.round(h*ratio);
 var ctx=cv.getContext("2d");ctx.setTransform(ratio,0,0,ratio,0,0);ctx.clearRect(0,0,w,h);ctx.fillStyle="#fff";ctx.fillRect(0,0,w,h);
 return{ctx:ctx,w:w,h:h};
}

function drawLine(cv,series,opt){
 opt=opt||{};
 var s=canvasSetup(cv,opt.height||430),ctx=s.ctx,w=s.w,h=s.h,p={l:62,r:18,t:16,b:36};
 var xmin=opt.xmin==null?0:opt.xmin,xmax=opt.xmax==null?series[0].x[series[0].x.length-1]:opt.xmax;
 var ymin=Infinity,ymax=-Infinity;
 series.forEach(function(a){
  if(a.hide)return;
  for(var i=0;i<a.y.length;i++){
   var xx=a.x[i];if(xx<xmin||xx>xmax)continue;
   var yy=a.y[i];if(Number.isFinite(yy)){ymin=Math.min(ymin,yy);ymax=Math.max(ymax,yy)}
  }
 });
 if(opt.ymin!=null)ymin=opt.ymin;if(opt.ymax!=null)ymax=opt.ymax;
 if(!Number.isFinite(ymin)||!Number.isFinite(ymax)){ymin=0;ymax=1}
 if(ymin===ymax){ymin-=1;ymax+=1}
 if(opt.ymin==null){var m=(ymax-ymin)*0.05||1;ymin-=m}
 if(opt.ymax==null){var m2=(ymax-ymin)*0.05||1;ymax+=m2}
 var X=function(x){return p.l+(x-xmin)/(xmax-xmin)*(w-p.l-p.r)};
 var Y=function(y){return p.t+(ymax-y)/(ymax-ymin)*(h-p.t-p.b)};
 ctx.font="11px system-ui";ctx.strokeStyle="#e5e7eb";ctx.fillStyle="#667085";ctx.lineWidth=1;
 for(var k=0;k<=5;k++){
  var yy=p.t+(h-p.t-p.b)*k/5;ctx.beginPath();ctx.moveTo(p.l,yy);ctx.lineTo(w-p.r,yy);ctx.stroke();
  ctx.fillText(fmt(ymax-(ymax-ymin)*k/5,opt.out?0:1),5,yy+4);
 }
 for(k=0;k<=6;k++){
  var xx=p.l+(w-p.l-p.r)*k/6;ctx.beginPath();ctx.moveTo(xx,p.t);ctx.lineTo(xx,h-p.b);ctx.stroke();
  ctx.fillText(fmt((xmin+(xmax-xmin)*k/6)/60,0),xx-8,h-11);
 }
 ctx.fillText("Time [min]",w/2-26,h-3);
 ctx.save();ctx.beginPath();ctx.rect(p.l,p.t,w-p.l-p.r,h-p.t-p.b);ctx.clip();
 series.forEach(function(a){
  if(a.hide)return;ctx.strokeStyle=a.color;ctx.globalAlpha=a.alpha==null?1:a.alpha;ctx.lineWidth=a.width||2;ctx.setLineDash(a.dash||[]);
  ctx.beginPath();var started=false,stride=Math.max(1,Math.floor(a.x.length/Math.max(1000,w*2)));
  for(var i=0;i<a.x.length;i+=stride){
   var xv=a.x[i];if(xv<xmin||xv>xmax)continue;var px=X(xv),py=Y(a.y[i]);
   if(!started){ctx.moveTo(px,py);started=true}else ctx.lineTo(px,py);
  }
  ctx.stroke();
 });
 ctx.restore();ctx.globalAlpha=1;ctx.setLineDash([]);ctx.strokeStyle="#98a2b3";ctx.strokeRect(p.l,p.t,w-p.l-p.r,h-p.t-p.b);
}

function drawDual(cv,temp,output,c){
 var s=canvasSetup(cv,450),ctx=s.ctx,w=s.w,h=s.h,p={l:62,r:54,t:16,b:36},xmin=2400,xmax=4200;
 var tmin=Infinity,tmax=-Infinity;
 for(var i=0;i<temp.x.length;i++)if(temp.x[i]>=xmin&&temp.x[i]<=xmax){tmin=Math.min(tmin,temp.y[i]);tmax=Math.max(tmax,temp.y[i])}
 var mg=(tmax-tmin)*0.12||10;tmin-=mg;tmax+=mg;
 var X=function(x){return p.l+(x-xmin)/(xmax-xmin)*(w-p.l-p.r)};
 var YT=function(y){return p.t+(tmax-y)/(tmax-tmin)*(h-p.t-p.b)};
 var YO=function(y){return p.t+(c.umax-y)/(c.umax-c.umin)*(h-p.t-p.b)};
 ctx.font="11px system-ui";ctx.strokeStyle="#e5e7eb";ctx.fillStyle="#667085";
 for(var k=0;k<=5;k++){
  var yy=p.t+(h-p.t-p.b)*k/5;ctx.beginPath();ctx.moveTo(p.l,yy);ctx.lineTo(w-p.r,yy);ctx.stroke();
  ctx.fillText(fmt(tmax-(tmax-tmin)*k/5,1),5,yy+4);
  ctx.fillText(fmt(c.umax-(c.umax-c.umin)*k/5,0)+"%",w-p.r+6,yy+4);
 }
 for(k=0;k<=6;k++){var xx=p.l+(w-p.l-p.r)*k/6;ctx.beginPath();ctx.moveTo(xx,p.t);ctx.lineTo(xx,h-p.b);ctx.stroke();ctx.fillText(fmt((xmin+(xmax-xmin)*k/6)/60,0),xx-8,h-11)}
 ctx.fillText("Time [min]",w/2-26,h-3);
 ctx.save();ctx.beginPath();ctx.rect(p.l,p.t,w-p.l-p.r,h-p.t-p.b);ctx.clip();
 ctx.strokeStyle=colors.onoff;ctx.lineWidth=2.2;ctx.beginPath();var started=false;
 for(i=0;i<temp.x.length;i++){if(temp.x[i]<xmin||temp.x[i]>xmax)continue;var px=X(temp.x[i]),py=YT(temp.y[i]);if(!started){ctx.moveTo(px,py);started=true}else ctx.lineTo(px,py)}ctx.stroke();
 ctx.strokeStyle=colors.output;ctx.lineWidth=1.8;ctx.beginPath();started=false;var lastY=0;
 for(i=0;i<output.x.length;i++){if(output.x[i]<xmin||output.x[i]>xmax)continue;px=X(output.x[i]);py=YO(output.y[i]);if(!started){ctx.moveTo(px,py);started=true}else{ctx.lineTo(px,lastY);ctx.lineTo(px,py)}lastY=py}ctx.stroke();
 ctx.restore();ctx.strokeStyle="#98a2b3";ctx.strokeRect(p.l,p.t,w-p.l-p.r,h-p.t-p.b);
}

function paramLabel(kind,p){
 if(kind==="p")return "Kp="+fmt(p.kp,3);
 if(kind==="pi")return "Kc="+fmt(p.kc,3)+", Ti="+fmt(p.ti,0)+"s";
 return "Kc="+fmt(p.kc,3)+", Ti="+fmt(p.ti,0)+"s, Td="+fmt(p.td,0)+"s";
}

function renderHistoryLegend(kind,current){
 var el=$(kind+"legend"),parts=[];
 parts.push('<span><i class="dot" style="background:'+colors.base+'"></i>기본 '+paramLabel(kind,fixed[kind])+'</span>');
 histories[kind].forEach(function(h,i){parts.push('<span><i class="dot" style="background:'+historyColors[i%historyColors.length]+'"></i>'+h.label+'</span>')});
 parts.push('<span><i class="dot" style="background:'+colors[kind]+'"></i>현재 '+paramLabel(kind,current)+'</span>');
 el.innerHTML=parts.join("");
}

function drawTuning(kind,c,currentResult,baseResult){
 var params=currentParams(kind),series=[{x:baseResult.t,y:baseResult.y,color:colors.base,dash:[7,5],width:1.5,alpha:.8}];
 histories[kind].forEach(function(h,i){series.push({x:h.r.t,y:h.r.y,color:historyColors[i%historyColors.length],width:1.5,alpha:.72})});
 series.push({x:currentResult.t,y:currentResult.y,color:colors[kind],width:2.7});
 drawLine($(kind+"chart"),series,{height:360});
 renderHistoryLegend(kind,params);
 var m=perf(c,currentResult);
 $(kind+"final").textContent=fmt(m.fin,2)+" ℃";
 $(kind+"err").textContent=fmt(m.err,2)+" ℃";
 $(kind+"over").textContent=fmt(m.ov,2)+" %";
}

function saveTrace(kind){
 try{
  var c=getConfig(),p=currentParams(kind),r=simulate(c,kind,p);
  histories[kind].push({label:paramLabel(kind,p),r:r});
  if(histories[kind].length>6)histories[kind].shift();
  runAll();
 }catch(e){showError(e)}
}

function clearHistories(){
 histories={p:[],pi:[],pid:[]};
}

function showError(e){
 var st=$("status");st.className="status err";st.textContent=e.message||String(e);
}

function runAll(){
 var st=$("status");
 try{
  st.className="status";st.textContent="계산 중…";
  var c=getConfig();
  var baseOn=simulate(c,"onoff",{});
  var baseP=simulate(c,"p",fixed.p);
  var basePI=simulate(c,"pi",fixed.pi);
  var basePID=simulate(c,"pid",fixed.pid);
  var curP=simulate(c,"p",currentParams("p"));
  var curPI=simulate(c,"pi",currentParams("pi"));
  var curPID=simulate(c,"pid",currentParams("pid"));
  var osc=lastWindowWidth(c,baseOn,1800),pm=perf(c,baseP),pim=perf(c,basePI),pidm=perf(c,basePID);
  $("osc").textContent=fmt(osc.width,2)+" ℃";
  $("oscsub").textContent="min "+fmt(osc.min,2)+" / max "+fmt(osc.max,2);
  $("pres").textContent=fmt(pm.err,2)+" ℃";$("pressub").textContent="최종 "+fmt(pm.fin,2)+" ℃";
  $("piderr").textContent=fmt(pidm.err,6)+" ℃";$("piderrsub").textContent="최종 "+fmt(pidm.fin,2)+" ℃";
  $("pierr").textContent=fmt(pim.err,4)+" ℃";$("pierrsub").textContent="추가 비교용";
  drawLine($("onoffFull"),[{x:baseOn.t,y:baseOn.y,color:colors.onoff,width:2.2}],{height:430});
  drawDual($("onoffZoom"),{x:baseOn.t,y:baseOn.y},{x:baseOn.t,y:baseOn.u},c);
  var spArr=new Float64Array(baseOn.t.length);spArr.fill(c.sp);
  drawLine($("compareChart"),[
   {x:baseOn.t,y:baseOn.y,color:colors.onoff,width:2},
   {x:baseP.t,y:baseP.y,color:colors.p,width:2},
   {x:basePID.t,y:basePID.y,color:colors.pid,width:2},
   {x:baseOn.t,y:spArr,color:colors.sp,dash:[7,5],width:1.4}
  ],{height:470,ymin:750,ymax:1100});
  drawTuning("p",c,curP,baseP);
  drawTuning("pi",c,curPI,basePI);
  drawTuning("pid",c,curPID,basePID);
  lastBundle={c:c,onoff:baseOn,p:curP,pi:curPI,pid:curPID};
  st.textContent="계산 완료 · 기본 과제 결과 + 현재 튜닝 결과가 갱신되었습니다.";
 }catch(e){showError(e)}
}

function schedule(){clearTimeout(timer);timer=setTimeout(runAll,55)}

function bindTune(rangeId,numId,kind){
 var r=$(rangeId),n=$(numId);
 r.addEventListener("input",function(){n.value=r.value;schedule()});
 n.addEventListener("input",function(){var v=Number(n.value);if(Number.isFinite(v))r.value=clamp(v,Number(r.min),Number(r.max));schedule()});
 r.addEventListener("change",function(){saveTrace(kind)});
 n.addEventListener("change",function(){saveTrace(kind)});
}

function bindProcess(id,eventName){
 $(id).addEventListener(eventName||"input",function(){
  clearHistories();
  if(id==="model"){var custom=this.value==="custom";$("fopdt").classList.toggle("hidden",custom);$("custom").classList.toggle("hidden",!custom)}
  schedule();
 });
}

bindTune("pk","pkn","p");
bindTune("pikc","pikcn","pi");bindTune("piti","pitin","pi");
bindTune("pidkc","pidkcn","pid");bindTune("pidti","pidtin","pid");bindTune("pidtd","pidtdn","pid");

["initial","bias","tau","gain","delay","cinitial","cdelay","equation","sp","duration","dt","umin","umax","hys"].forEach(function(id){bindProcess(id,"input")});
bindProcess("model","change");bindProcess("aw","change");

$("clearHistory").addEventListener("click",function(){clearHistories();runAll()});
$("reset").addEventListener("click",function(){
 Object.keys(defaults).forEach(function(k){if($(k))val(k,defaults[k])});
 val("model","fopdt");val("aw","on");
 val("pk",defaults.pk);val("pkn",defaults.pk);
 val("pikc",defaults.pikc);val("pikcn",defaults.pikc);val("piti",defaults.piti);val("pitin",defaults.piti);
 val("pidkc",defaults.pidkc);val("pidkcn",defaults.pidkc);val("pidti",defaults.pidti);val("pidtin",defaults.pidti);val("pidtd",defaults.pidtd);val("pidtdn",defaults.pidtd);
 $("fopdt").classList.remove("hidden");$("custom").classList.add("hidden");clearHistories();runAll();
});
$("csv").addEventListener("click",function(){
 if(!lastBundle)return;
 var b=lastBundle,s="time_s,onoff_T,onoff_u,P_T,P_u,PI_T,PI_u,PID_T,PID_u\n";
 for(var i=0;i<b.onoff.t.length;i++)s+=b.onoff.t[i]+","+b.onoff.y[i]+","+b.onoff.u[i]+","+b.p.y[i]+","+b.p.u[i]+","+b.pi.y[i]+","+b.pi.u[i]+","+b.pid.y[i]+","+b.pid.u[i]+"\n";
 var url=URL.createObjectURL(new Blob([s],{type:"text/csv"})),a=document.createElement("a");a.href=url;a.download="control_simulation.csv";a.click();setTimeout(function(){URL.revokeObjectURL(url)},500);
});
window.addEventListener("resize",schedule);
runAll();
})();
