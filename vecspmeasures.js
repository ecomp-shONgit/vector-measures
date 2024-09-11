/**

    2022 Prof. Schäfer Uni Trier, Prof. Schubert Uni Leipzig
    Measures on vector represented mathematical objects (to compute distances between objects in a vector space).

**/
"use strict";

/* !equal size of vector v1 and v2! */

/* pythangoraeic difference */
function euclideanM( v1, v2 ){
    let d = 0.0;
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    for(let i = 0; i < l; i += 1 ){
        let z = v1[i] - v2[i];
        d += z*z;
    }
    return Math.sqrt( d );

}

function L2( x1, x2 ){
    return euclideanM( x1, x2 );
}

function chebyshevM( v1, v2 ){
    let da = [];
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    for(let i = 0; i < l; i += 1 ){
        da.push( Math.abs( v1[i] - v2[i] ) );
    }
    return  Math.max( ...da );
}

function chebyaddinterseczM( v1, v2 ){
    return chebyshevM( v1, v2 )+(intersectionZMN(v1, v2)*255); //250 ????
}

function minkowskiM( v1, v2, order ){
    let d = 0.0;
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    if( order < 0.0 ){
        return NaN;
    }
    for(let i = 0; i < l; i += 1 ){
        d += Math.pow( Math.abs( v1[i] - v2[i] ), order );
    }
    return Math.pow( d, (1.0/order) );
}

function manhattanM( v1, v2 ){ //also city block distance, rectilinear distance, taxicab norm
    let d = 0.0;
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    for(let i = 0; i < l; i += 1 ){
        d +=  Math.abs( v1[i] - v2[i] );
    }
    return d;
}



/*absolut difference*/
function canberraM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for(let i = 0; i < l; i += 1 ){
        let temp = Math.abs(v1[i])  + Math.abs(v2[i]);
        if(temp != 0 && temp != 0.0){
            d += Math.abs(v1[i]  - v2[i]) / temp;
        }
    }
    
    return d;
}

function soerensenM( v1, v2 ){ //also Bray-Curtis
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let a = 0.0;
    let b = 0.0
    for( let i = 0; i < l; i += 1 ){
        a += Math.abs(v1[i]  - v2[i]);
        b += Math.abs(v1[i])  + Math.abs(v2[i]); //is this abs right?
    }
    return a/b;
}

function normvecM( v1, v2 ){ 
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
	//vektoren normieren
	let n1 = 0.0;
	let n2 = 0.0;
	for( let i = 0; i < l; i += 1 ){
	    n1 += v1[i]*v1[i];
	
	    n2 = v2[i]*v2[i];
    }
    n1 = Math.sqrt( n1 );
    n2 = Math.sqrt( n2 );


	
	//euklidian kombination
	let normval = 0.0
    for( let i = 0; i < l; i += 1 ){
        let tempnor1 = v1[i]/n1;
        let tempnor2 = v2[i]/n2;
        let tempdiff = (tempnor1 - tempnor2);
        normval += (tempdiff*tempdiff);
    }
    normval = normval / l;
    let retval = Math.sqrt(normval)  || 0.0;
    
    return parseFloat( retval);
    /*let normval = 0.0
    for( let i = 0; i < l; i += 1 ){
        let tempnor1 = v1[i]/n1;
        let tempnor2 = v2[i]/n2;
        let tempdiff = tempnor1 - tempnor2;
        normval += Math.abs(tempdiff*tempdiff);
    }
    normval = normval / l;
    
    return parseFloat( normval);
    */
}

function gowerM( v1, v2, range ){ //range of values in the row, i.e. same feature; i.e. special case all features have the same potential range
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    let r = 0.0;
    for( let i = 0; i < l; i += 1 ){
        let ran = 1.0;
        if( v1[i] == 0 || v2[i] == 0){ 
            ran = 0.0; //general absens or presens
        }
        //d += (Math.abs(v1[i]  - v2[i])/Math.max(v1[i], v2[i]) ); //unclude a scaling for the range of both values
        d += (Math.abs(v1[i]  - v2[i])/range );
        r += ran;
    }
    return (d/r);
}

function soergelM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let a = 0.0;
    let b = 0.0;
    for( let i = 0; i < l; i += 1 ){
        a += Math.abs(v1[i]  - v2[i]);
        b += Math.max(v1[i], v2[i]); 
    }
    return a/b;
}

/*function kulczynskiM( v1, v2 ){ //binary purpose
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let a = 0.0;
    let b = 0.0;
    let c = 0.0
    for( let i = 0; i < l; i += 1 ){
        //a += Math.abs(v1[i]  - v2[i]);
        //b += Math.min(v1[i], v2[i]); 
        a += Math.min(v1[i], v2[i]); 
        b += v1[i];
        c += v2[i];
    }
    //return a/b;
    return 0.5*((a/b)+(a/c));
}*/

function lorentzianM( v1, v2 ){ 
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        d += Math.log( 1 + Math.abs(v1[i]  - v2[i]) );
    }
    return d;
}
/*intersections*/

function intersectionM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        d += Math.min(v1[i], v2[i]);
    }
    return d;
}
function intersectionZM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    let sv1 = 0.0;
    let sv2 = 0.0;
    for( let i = 0; i < l; i += 1 ){
        d += Math.min(v1[i], v2[i]);
        sv1 += v1[i];
        sv2 += v2[i];
    }
     
    return (1-(d/Math.min(sv1, sv2)));
}

function intersectionZMN( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 1.0;
    let sv1 = 1.0;
    let sv2 = 1.0;
    let D = 1.0;
    for( let i = 0; i < l; i += 1 ){
        d += Math.min(v1[i], v2[i]);
        sv1 += v1[i];
        sv2 += v2[i];
        D += Math.max(v1[i], v2[i]);
    }
    if(D == 0.0){
        D = 1.0;
    }
    let minofsv = Math.min(sv1, sv2);
    if(minofsv == 0.0 ){
        minofsv = 1;
    }
    return (1-(d/minofsv))*D*0.001;
}

function wavehedgesM( v1, v2 ){ //wave edges; Edge waves from fluid dynamics ??? - introduced by Chan 2007 very bad, no sources but many modern replications
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        //console.log(d, Math.abs(v1[i]  - v2[i]), Math.max(v1[i], v2[i]), Math.abs(v1[i]  - v2[i])/Math.max(v1[i], v2[i]));
        const u = Math.max( v1[i], v2[i] );
        if( 0 != u ){ //get zero values not wrong 
            d += Math.abs(v1[i]  - v2[i]) / u; 
        }
        //changed acording to: https://arxiv.org/pdf/1409.0923.pdf Hassant
        //d += 1-( ( 1 + mi ) / ( 1 + ma ) );
    }
    //console.log(d);
    return d;
}

function hassanatM( v1, v2 ){ //https://arxiv.org/pdf/1409.0923.pdf
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        let mi = Math.min(v1[i], v2[i]);
        let ma = Math.max(v1[i], v2[i]);
        let diffab = Math.abs( v1[i] - v2[i]);
        if( mi >= 0 ){
            d +=  ( ( diffab ) / ( 1 + ma ) );
        } else {
            
            let mia = Math.abs( mi );
            d += ( ( diffab ) / ( 1 + ma + mia ) );    
        }
    }
    return d;
}

//Sobolev https://www.techfak.uni-bielefeld.de/~fschleif/mlr/mlr_03_2007.pdf

function motykaM( v1, v2 ){ 
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let a = 0.0;
    let b = 0.0;
    for( let i = 0; i < l; i += 1 ){
        a += v1[i]  + v2[i];
        b += Math.max(v1[i], v2[i]); 
    }
    return b/a;
}

/*function kulczynskisM( v1, v2 ){ // see above and binary data 
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let a = 0.0;
    let b = 0.0;
    for( let i = 0; i < l; i += 1 ){
        a += Math.abs(v1[i]  - v2[i]);
        b += Math.min(v1[i], v2[i]); 
    }
    return (a/b);
}*/

function ruzickaM( v1, v2 ){ // 1- tanimotoM; minmax
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let a = 0.0;
    let b = 0.0;
    for( let i = 0; i < l; i += 1 ){
        a += Math.min(v1[i], v2[i]);
        b += Math.max(v1[i], v2[i]); 
    }
    return a/b;
}

function tanimotoM( v1, v2 ){ //also Jaccard
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let a = 0.0;
    let b = 0.0;
    for( let i = 0; i < l; i += 1 ){
        let ma = Math.max(v1[i], v2[i]);
        a += ma - Math.min(v1[i], v2[i]); 
        b += ma;
    }
    return a/b;
}


/*inner product*/
function innerproductM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        d += v1[i] * v2[i];
    }
    return d;
    
}

function harmonicmeanM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        const u = (v1[i] + v2[i])
        if( u != 0 ){
            d += (v1[i] * v2[i]) / u;
        }
    }
    return d*2;
    
}

function cosineM( v1, v2 ){ 
    
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    
    let n1 = 0.0;
    let n2 = 0.0;
    let pu = 0.0;
    for(let i = 0; i < l; i += 1 ){
        n1 += (v1[i]+1)*(v1[i]+1);
        n2 += (v2[i]+1)*(v2[i]+1);
        pu += (v1[i]+1)*(v2[i]+1);
    }
    n1 = Math.sqrt(n1);
    n2 = Math.sqrt(n2);
    let re =  1.0 - (pu / (n1*n2)); 
    if(!re){ // that means therer is a NaN input or a zero vector
        //console.log("NaN", v1, v2);
        re = 0.0;
    }
    return  re;
}

function kumarhassebrookM( v1, v2 ){ //peak to correlation energy PCE
    
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let n1 = 0.0;
    let n2 = 0.0;
    let pu = 0.0;
    for(let i = 0; i < l; i += 1 ){
        n1 += v1[i]*v1[i];
        n2 += v2[i]*v2[i];
        pu += v1[i]*v2[i];
    }
    
    let re =  1.0 - (pu / (n1+n2-pu)); 
    if(!re){ // that means therer is a NaN input 
        re = 0.0;
    }
    return  re;
}

function diceM( v1, v2 ){ 
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let n1 = 0.0;
    let n2 = 0.0;
    let pu = 0.0;
    for(let i = 0; i < l; i += 1 ){
        n1 += v1[i]*v1[i];
        n2 += v2[i]*v2[i];
        pu += v1[i]*v2[i];
    }
    
    let re =  1.0 - (pu / (n1+n2)); 
    if(!re){ // that means therer is a NaN input 
        re = 0.0;
    }
    return  re;
}

/*fiedelity based*/
function fidelityM( v1, v2 ){ //sum of geometric means
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for(let i = 0; i < l; i += 1 ){
        d += Math.sqrt( v1[i]*v2[i] );
        
    }
    return  d;
}

function bhattacarya1M( v1, v2 ){ 
    const zwierg = fidelityM( v1, v2 );
    let z = 0;
    if( zwierg > 1.0 ){
        z = Math.acos( 1/zwierg );//for values greater than 1
    } else {
        z = Math.acos( zwierg );
    }
    return  z*z;
}

function bhattacarya2M( v1, v2 ){
    //return  -1*Math.log( fidelityM( v1, v2 ) );
    return  Math.abs( Math.log( fidelityM( v1, v2 ) ) ); //for values greater that 1
}

function hellingerM( v1, v2 ){
    //return 2*Math.sqrt( 1-fidelityM( v1, v2 ) );
    return 2*Math.sqrt( fidelityM( v1, v2 ) ); //for values greater than 1 e.g. not frequnecies
}


/*chisquared variants*/

function squaredeuclideanM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for(let i = 0; i < l; i += 1 ){
        let z = v1[i]-v2[i];
        d += z * z;
        
    }
    return  d;
}

function pearsonchisquaredM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for(let i = 0; i < l; i += 1 ){
        if(v2[i] != 0){
            let z = v1[i]-v2[i];
            
            d += (z * z) / v2[i];
        }
    }
    return  d;
}

function neymanchisquaredM( v1, v2 ){ 
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for(let i = 0; i < l; i += 1 ){
        let z = v1[i]-v2[i];
        if( z != 0 && v1[i] != 0 ){
            d += (z * z) / v1[i];
        }
        
    }
    return  d;
}

function squaredchisquaredM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for(let i = 0; i < l; i += 1 ){
        let z = v1[i] - v2[i];
        let y = (v1[i]+v2[i]);
        if( z != 0 && y != 0 ){
            d += (z * z) / y;
        }
        
    }
    return  d;
}

function divergenceM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for(let i = 0; i < l; i += 1 ){
        let z = 1+(v1[i] - v2[i]);//if part of the vector is zero
        let y = 1+(v1[i] + v2[i]);
        d += (z * z) / (y * y); 
        
    }
    return  d;
}

function clarckM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        if( v1[i] != 0 && v2[i] != 0 ){
            let z = Math.abs(v1[i] - v2[i]) / (v1[i] + v2[i]);
            d += z * z ;
        }
    }
    return  Math.sqrt( d );
}

function additivesymmetricchisquaredM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        if( v1[i] != 0 && v2[i] != 0 ){
            let z = v1[i] - v2[i];
            let y = v1[i]+v2[i];
            d += ( (z * z) * y ) / (v1[i] * v2[i]);  
        }     
    }
    return  d;
}

/*shannon entropy alike*/
function kullbackleiblerM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        if( v1[i] != 0 && v2[i] != 0 ){
            d += v1[i] * Math.log( v1[i] / v2[i] ); 
        }
    }
    return  d;
}

function JeffreysM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        if( v1[i] != 0 && v2[i] != 0 ){
            d += ((v1[i] - v2[i]) * Math.log( v1[i] / v2[i] )  );
        } 
    }
    
    return  d;
}

function kullbackdivergenceM( v1, v2 ){ 
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        if( v1[i] != 0 && v2[i] != 0 ){
            
            d += v1[i] * Math.log( (2*v1[i]) / (v1[i]+v2[i]) ); 
        }
    }
    return  d;
}

function topsoeeM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        if( v1[i] != 0 && v2[i] != 0 ){
            d += (v1[i] * Math.log( (2*v1[i]) / (v1[i]+v2[i]) )) + (v2[i] * Math.log( (2*v2[i]) / (v1[i]+v2[i]) )); 
        }
    }
    return  d;
}

function jensenshannonM( v1, v2 ){ 
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let a = 0.0;
    let b = 0.0;
    for( let i = 0; i < l; i += 1 ){
        if( v1[i] != 0.0 && v2[i] != 0.0 ){
            a += (v1[i] * Math.log( (2*v1[i]) / (v1[i]+v2[i]) )); //(v1[i]+v2[i]) - mixture maybe expressed as: (v1[i]+v2[i])/2
            b += (v2[i] * Math.log( (2*v2[i]) / (v1[i]+v2[i]) )); 
        }
    }
    return  (0.5*(a+b) )  || 0.0;
}

function jensenM( v1, v2 ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        if( v1[i] != 0.0 && v2[i] != 0.0 ){
        let z = ((v1[i]+v2[i])/2);
            //console.log(z, v1[i], Math.log(v1[i]), v2[i],  Math.log(v2[i]), z*Math.log(z))
            d += (((( v1[i] * Math.log(v1[i]) ) + ( v2[i] * Math.log(v2[i]) ) ) / 2 ) - (z*Math.log(z)) );
        //console.log(d);
        } 
    }
    return  d  || 0.0;
}

/* stylo specific measures*/ 
function stdabw( VV ){ 
    let ll = VV.length; //anzahl der vectoren
    let l = VV[0].length;//anzahl der wörter pro vector
    let means = [];
    let stdabws = [];
    for( let i = 0; i < l; i+= 1){ //initalisieren mit der anzahl der wörter
        means.push( 0.0 );
        stdabws.push( 0.0 );
    }
    for( let i = 0; i < l; i += 1 ){//ech word form
        for( let j = 0; j < ll; j += 1 ){//each frequency value
            means[i] += VV[j][i];
        }
        means[i] = means[i] / ll;
    }

    for( let i = 0; i < l; i += 1 ){//ech word form
        for( let j = 0; j < ll; j += 1 ){//each frequency value
            stdabws[i] += ( means[i] - VV[j][i] ) * ( means[i] - VV[j][i] );
        }
        stdabws[i] = Math.sqrt( stdabws[i] / ll );
    }
    return [stdabws, means];
}


function edersimpleM( v1, v2 ){
    
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let v11 = [];
    let v22 = [];
    for( let i = 0; i < l; i += 1 ){
        v11.push( Math.sqrt( v1[i] ) );
        v22.push( Math.sqrt( v2[i] ) );
    }
    return manhattanM( v11, v22 );
}

function burrowsdeltaM( v1, v2, stami ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    //console.log(sta);
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        //d += Math.abs(  v1[i] - v2[i] ) / (stami[1][i]+1) ; 
        d += Math.abs( ( ( v1[i] - stami[1][i] ) / ( stami[0][i] + 1 ) ) - ( ( v2[i] - stami[1][i] ) / ( stami[0][i] + 1 ) ) );
    }
    //d /= l;
    return d;
}

function burrowsdeltaFastM( v1, v2, stami ){
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    //console.log(sta);
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        d += Math.abs(  v1[i] - v2[i] ) / (stami[0][i]+1) ; 
        //d += Math.abs( ( ( v1[i] - stami[1][i] ) / ( stami[0][i] + 1 ) ) - ( ( v2[i] - stami[1][i] ) / ( stami[0][i] + 1 ) ) );
    }
    d /= l;
    return d;
}

function argamonlineardeltaM( v1, v2, stami ){ //check definition 
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        let t = v1[i] - v2[i];
        d += Math.sqrt( ( t * t ) / (stami[1][i]+1) ); 
        
    }
    d /= l;
    return d;
}

function edersdeltaM( v1, v2, stami ){ 
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        const rankscalar = ( ( ( l - i ) + 1 ) / l );
        //d +=  Math.abs( ( v1[i] - v2[i] ) / (stami[i]+1) ) * ( ( ( l - i ) + 1 ) / l ); 
        d += Math.abs( ( ( v1[i] - stami[1][i] ) / ( stami[0][i] + 1 ) * rankscalar ) - ( ( v2[i] - stami[1][i] ) / ( stami[0][i] + 1 ) * rankscalar ) );
    }
    //d /= l;
    return d;
}

function edersdeltaFastM( v1, v2, stami ){ 
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        //const rankscalar = ( ( ( l - i ) + 1 ) / l );
        d +=  Math.abs( ( v1[i] - v2[i] ) / (stami[0][i]+1) ) * ( ( ( l - i ) + 1 ) / l ); 
        //d += Math.abs( ( ( v1[i] - stami[1][i] ) / ( stami[0][i] + 1 ) * rankscalar ) - ( ( v2[i] - stami[1][i] ) / ( stami[0][i] + 1 ) * rankscalar ) );
    }
    d /= l;
    return d;
}

function argamonsquadraticdeltaM( v1, v2, stami ){ 
    let l = v1.length; 
    if( l != v2.length ){
        return NaN;
    }
    let d = 0.0;
    for( let i = 0; i < l; i += 1 ){
        let t = v1[i] - v2[i];
        d += ( t * t ) / (stami[0][i]+1) ; 
    }
    d /= l;
    return d;
}

/*Transporttttt*/
function wassersorded(d, h, magni){
    let LL = Math.min( d.length, h.length ); //need to be cared for lengthyst array
    let wheretoput = 0;
    let i = 0;
    let j = 0;
    let tocarrayatthemoment = d[i];
    let platz = h[j];
    let carryon = true;

    let flow = 1;
    let work = 1;


    while( carryon ){
        
        //compute
        platz = (tocarrayatthemoment - platz) * (-1);
        //console.log(platz, i, j);
        
        if( platz > 0.0 ){ //can put all things in and have more space
            
            flow += tocarrayatthemoment;
            work += tocarrayatthemoment * ( Math.abs( i-j ) );
        } else if( platz < 0.0 ){ //cant put all in have more work to do 
            flow += d[i]-platz;
            work += ( d[i] - platz ) * ( Math.abs( i-j ) );
        } else {
            flow += tocarrayatthemoment;
            work += tocarrayatthemoment * ( Math.abs( i-j ) );
        }

        //check
        i += 1;
        if( !(i < LL) ){
            carryon = false;
            continue;
        }

        //advanc
        if( platz > 0.0 ){ //can put all things in and have more space
            tocarrayatthemoment = d[i];
        } else if( platz < 0.0 ){ //cant put all in have more work to do 
            tocarrayatthemoment = d[i] - (platz*1);
            j += 1;
            platz = h[j];
        } else {
            tocarrayatthemoment = d[i];
            j += 1;
            platz = h[j];
        }
        //console.log(platz, i, j);
    }

    
    let w1dist = work / flow;
    //rest
    platz = Math.abs(platz);
    if( platz != 0 ){
        w1dist += (platz* magni); //remainning dirt / place frome the last
    }
    //compute constant offset of uncompared parts
    if( j < h.length ){
        let add = 0;
        for( let t = j; t < h.length; t += 1 ){
            add += h[t]*magni;
        }
        w1dist += ( add / (Math.abs(j-h.length) ) );
    }
    if( i < i.length ){
        let add = 0;
        for( let t = i; t < d.length; t += 1 ){
            add += d[t]*magni;
        }
        w1dist += ( add / (Math.abs(i-d.length) ) );
    }

    return w1dist;

}

function wasserst1dM( d, h, magni ){
    //seyetrical comarison
    return ((wassersorded(d, h, magni)+wassersorded(h, d, magni))/2);
    
}

//
function skinhorn(){

}

/*TESTSW*/
function testmeasagainstAmeasure( m1, ms, ms2, a1, bs ){
    let resultingdiffs = [];
    let diffs1 = [];
    let stdif = stdabw( bs );
    //console.log("Standadabweichung:", stdif);
    for( let b = 0; b < bs.length; b += 1 ){
        const fktname = m1.name; 
        if( fktname != "minkowskiM" && 
            fktname != "burrowsdeltaM" && 
            fktname != "argamonlineardeltaM" && 
            fktname != "edersdeltaM" && 
            fktname != "argamonsquadraticdeltaM" &&
            fktname != "wasserst1dM" &&
            fktname != "gowerM" &&
            fktname != "edersdeltaFastM" &&
            fktname != "burrowsdeltaFastM"){
                diffs1.push( m1(a1, bs[b]) );
            } else {
                diffs1.push( m1(a1, bs[b], stdif ) );
            }
    }
    for( let m = 0; m < ms.length; m += 1 ){
        let rrr = [];
        for( let b = 0; b < bs.length; b += 1 ){
            //const er =  Math.abs(diffs1[b] - ms[m]( a1,bs[b] ) );
            //const er = (ms[m]( a1,bs[b] ) / diffs1[b]);
            
            const er = Math.abs(diffs1[b] - ms[m]( a1,bs[b] ) );
            rrr.push( er );
        }
        resultingdiffs.push( rrr );
    }
    for( let m = 0; m < ms2.length; m += 1 ){
        let rrr = [];
        for( let b = 0; b < bs.length; b += 1 ){
            //const er = Math.abs(diffs1[b] - ms2[m]( a1,bs[b], stdif ));
            //const er = (ms2[m]( a1,bs[b], stdif )  / diffs1[b]);
            const er = Math.abs( diffs1[b] - ms2[m]( a1,bs[b], stdif ) );
            /*if(!er){
                console.log( stdif, ms2[m]( a1,bs[b], stdif ), diffs1[b] );
            }*/
            rrr.push( er );
        }
        resultingdiffs.push( rrr );
    }
    return resultingdiffs;
}

function testvecmesure(){
    let A = [0.4, 54.9, 4.3];
    let B = [10.0, 1.9, 99.2];
    console.log("euclideanM: ", euclideanM(A,B));
    console.log("chebyshevM: ", chebyshevM(A,B));
    console.log("minkowskiM 0.5: ", minkowskiM(A,B, 0.5));
    console.log("manhattanM: ", manhattanM(A,B));
    

    console.log("canberraM: ", canberraM(A,B));
    console.log("soerensenM: ", soerensenM(A,B));
    console.log("normvecM: ", normvecM( A, B ), normvecM( B, A ));
    console.log("gowerM: ", gowerM(A,B, 100.0));
    console.log("soergelM: ", soergelM(A,B));
    //console.log("kulczynskiM: ", kulczynskiM(A,B));
    console.log("lorentzianM: ", lorentzianM(A,B));

    console.log("intersectionM: ", intersectionM(A,B));
    console.log("wavehedgesM: ", wavehedgesM(A,B));
    console.log("hassanatM: ",  hassanatM(A,B));

    console.log("motykaM: ", motykaM(A,B));
    //console.log("kulczynskisM: ", kulczynskisM(A,B));
    console.log("ruzickaM: ", ruzickaM(A,B));
    console.log("tanimotoM: ", tanimotoM(A,B));

    console.log("innerproductM: ", innerproductM(A,B));
    console.log("harmonicmeanM: ", harmonicmeanM(A,B));
    console.log("cosineM: ", cosineM(A,B));
    console.log("kumarhassebrookM: ", kumarhassebrookM(A,B));
    console.log("diceM: ", diceM(A,B));

    console.log("fidelityM: ", fidelityM(A,B));
    console.log("bhattacarya1M: ", bhattacarya1M(A,B));
    console.log("bhattacarya2M: ", bhattacarya2M(A,B));
    console.log("hellingerM: ",    hellingerM(A,B));

    console.log("jensenM: ", jensenM(A,B));
    console.log("jensenshannonM: ", jensenshannonM(A,B));
    console.log("topsoeeM: ", topsoeeM(A,B));
    console.log("kullbackdivergenceM: ",    kullbackdivergenceM(A,B));
    console.log("JeffreysM: ", JeffreysM(A,B));
    console.log("kullbackleiblerM: ", kullbackleiblerM(A,B));

    console.log("squaredeuclideanM: ", squaredeuclideanM(A,B));
    console.log("pearsonchisquaredM: ", pearsonchisquaredM(A,B));
    console.log("neymanchisquaredM: ", neymanchisquaredM(A,B));
    console.log("squaredchisquaredM: ",    squaredchisquaredM(A,B));
    console.log("divergenceM: ", divergenceM(A,B));
    console.log("clarckM: ", clarckM(A,B));      
    console.log("additivesymmetricchisquaredM: ", additivesymmetricchisquaredM(A,B)); 
           
    console.log("edersimpleM: ", edersimpleM(A,B));
    console.log("burrowsdeltaM: ",    burrowsdeltaM(A,B, stdabw([A,B])));
    console.log("argamonlineardeltaM: ", argamonlineardeltaM(A,B, stdabw([A,B])));
    console.log("edersdeltaM: ", edersdeltaM(A,B, stdabw([A,B])));      
    console.log("argamonsquadraticdeltaM: ", argamonsquadraticdeltaM(A,B,stdabw([A,B]))); 

    console.log("wasserst1dM", wasserst1dM(A,B, 0.01));
    console.log("wasserst1dM swi", wasserst1dM(B,A, 0.01));

    
}

//testvecmesure();


const dmeasuredict = {};
dmeasuredict["euclideanM"] = euclideanM;
dmeasuredict["chebyshevM"] = chebyshevM;
dmeasuredict["minkowskiM"] =  minkowskiM;
dmeasuredict["manhattanM"] =  manhattanM;
dmeasuredict["canberraM"] = canberraM;
dmeasuredict["soerensenM"] = soerensenM;
dmeasuredict["normvecM"] = normvecM;
dmeasuredict["gowerM"] = gowerM;
dmeasuredict["soergelM"] = soergelM;
dmeasuredict["lorentzianM"] = lorentzianM;
dmeasuredict["intersectionM"] = intersectionM;
dmeasuredict["wavehedgesM"] = wavehedgesM;
dmeasuredict["hassanatM"] = hassanatM;
dmeasuredict["motykaM"] = motykaM;
dmeasuredict["ruzickaM"] = ruzickaM;
dmeasuredict["tanimotoM"] = tanimotoM;
dmeasuredict["innerproductM"] = innerproductM;
dmeasuredict["harmonicmeanM"] = harmonicmeanM;
dmeasuredict["cosineM"] = cosineM;
dmeasuredict["kumarhassebrookM"] = kumarhassebrookM;
dmeasuredict["diceM"] = diceM;
dmeasuredict["fidelityM"] = fidelityM;
dmeasuredict["bhattacarya1M"] = bhattacarya1M;
dmeasuredict["bhattacarya2M"] = bhattacarya2M;
dmeasuredict["hellingerM"] = hellingerM;
dmeasuredict["jensenM"] = jensenM;
dmeasuredict["jensenshannonM"] = jensenshannonM;
dmeasuredict["topsoeeM"] = topsoeeM;
dmeasuredict["kullbackdivergenceM"] = kullbackdivergenceM;
dmeasuredict["JeffreysM"] = JeffreysM;
dmeasuredict["kullbackleiblerM"] = kullbackleiblerM;
dmeasuredict["squaredeuclideanM"] = squaredeuclideanM;
dmeasuredict["pearsonchisquaredM"] = pearsonchisquaredM;
dmeasuredict["neymanchisquaredM"] = neymanchisquaredM;
dmeasuredict["squaredchisquaredM"] = squaredchisquaredM;
dmeasuredict["divergenceM"] = divergenceM;
dmeasuredict["clarckM"] = clarckM;      
dmeasuredict["additivesymmetricchisquaredM"] = additivesymmetricchisquaredM; 
dmeasuredict["edersimpleM"] = edersimpleM;
dmeasuredict["burrowsdeltaM"] = burrowsdeltaM;
dmeasuredict["argamonlineardeltaM"] = argamonlineardeltaM;
dmeasuredict["edersdeltaM"] = edersdeltaM;      
dmeasuredict["argamonsquadraticdeltaM"] = argamonsquadraticdeltaM; 
dmeasuredict["wasserst1dM"] = wasserst1dM;

