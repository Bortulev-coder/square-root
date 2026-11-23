(function(root){
 'use strict';
 const DecimalLib = typeof module==='object' && module.exports ? require('./vendor/decimal.js') : root.Decimal;
 const D=DecimalLib.clone({precision:180,rounding:DecimalLib.ROUND_HALF_EVEN});
 function run({method='heron',a='2',x0='1',steps=12,tolerance='1e-12'}={}){
  if(!['heron','bakhshali','two'].includes(method))throw Error('Выберите один из трёх методов.');
  let rad,initial,tol;
  try{rad=new D(String(a).replace(',','.'));initial=method==='two'?rad:new D(String(x0).replace(',','.'));tol=new D(tolerance);}catch{throw Error('Введите числа: например 2, 0.25 или 1e-3.');}
  if(!rad.isFinite()||!rad.gt(0)||rad.lt('1e-50')||rad.gt('1e50'))throw Error('Для эксперимента введите a от 10⁻⁵⁰ до 10⁵⁰. В теории первые два метода допускают любое a > 0.');
  if(method==='two'&&!rad.lt(3))throw Error('Для метода с двумя переменными требуется 0 < a < 3. При a = 3 метод останавливается в нуле.');
  if(method!=='two'&&(!initial.isFinite()||!initial.gt(0)||initial.lt('1e-50')||initial.gt('1e50')))throw Error('Начальное приближение x₀ должно быть положительным: от 10⁻⁵⁰ до 10⁵⁰.');
  if(!Number.isInteger(steps)||steps<1||steps>60)throw Error('Число итераций должно быть целым от 1 до 60.');
  if(!tol.isFinite()||!tol.gt(0)||!tol.lt(1))throw Error('Порог ошибки должен находиться между 0 и 1.');
  const s=rad.sqrt(),rows=[];
  let x=method==='two'?rad:initial,c=rad.minus(1),first=null;
  for(let n=0;n<=steps;n++){
   const signed=x.div(s).minus(1),error=signed.abs(),residual=x.times(x).div(rad).minus(1).abs();
   if(first===null&&error.lte(tol))first=n;
   const row={n,x:x.toString(),error:error.toString(),signed:signed.toString(),residual:residual.toString(),c:method==='two'?c.toString():null,logError:error.isZero()?null:error.logarithm(10).toNumber()};
   rows.push(row);
   if(n===steps)break;
   if(method==='heron')x=x.plus(rad.div(x)).div(2);
   else if(method==='bakhshali'){
    const alpha=rad.minus(x.times(x)).div(x.times(2)),b=x.plus(alpha);
    row.alpha=alpha.toString();row.b=b.toString();
    x=b.minus(alpha.times(alpha).div(b.times(2)));
   }else{
    const nextX=x.times(new D(1).minus(c.div(2)));
    const nextC=c.times(c).times(c.minus(3)).div(4);
    x=nextX;c=nextC;
   }
  }
  return {method,a:rad.toString(),x0:rows[0].x,steps,tolerance:tol.toString(),root:s.toString(),rows,first,precision:180};
 }
 function format(s,digits=13){const d=new D(s);if(d.isZero())return '0';return d.toSignificantDigits(digits).toString();}
 const api={run,format,D};
 if(typeof module==='object'&&module.exports)module.exports=api;else root.SqrtEngine=api;
})(typeof globalThis!=='undefined'?globalThis:this);
