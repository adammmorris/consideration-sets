/* Statistics */

function getNormal(mean, std_dev) {
	var u = 2*Math.random()-1;
	var v = 2*Math.random()-1;
	var r = u*u + v*v;
	/*if outside interval [0,1] start over*/
	if(r == 0 || r > 1) return gaussRandom();

	var c = Math.sqrt(-2*Math.log(r)/r);
	return mean + u*c*std_dev;
}

function getRandomInt(min, max) {
    return Math.floor(Math.random() * (max - min + 1)) + min;
}

// spec should be: {value:weight, value:weight, ...}
function getRandomIntWeighted(spec) {
  var i, sum=0, r=Math.random();
  for (i in spec) {
    sum += spec[i];
    if (r <= sum) return i;
  }
}

/* Other */

function splice_better(array, index) {
    var x = new Array(index.length);
    for (var j = 0; j < index.length; j++) {
        x[j] = array[index[j]];
    }
    return x;
}

function shuffle(o){
    for(var j, x, i = o.length; i; j = Math.floor(Math.random() * i), x = o[--i], o[i] = o[j], o[j] = x);
    return o;
}

function save_data(data, table_name){
	$.ajax({
		type:'post',
		cache: false,
		url: 'savedata.php', // change this to point to your php file.
		// opt_data is to add additional values to every row, like a subject ID
		// replace 'key' with the column name, and 'value' with the value.
		data: {
			table: table_name,
			json: JSON.stringify(data),
		}
	});
}

function round(value, exp) {
  if (typeof exp === 'undefined' || +exp === 0)
    return Math.round(value);

  value = +value;
  exp = +exp;

  if (isNaN(value) || !(typeof exp === 'number' && exp % 1 === 0))
    return NaN;

  // Shift
  value = value.toString().split('e');
  value = Math.round(+(value[0] + 'e' + (value[1] ? (+value[1] + exp) : exp)));

  // Shift back
  value = value.toString().split('e');
  return +(value[0] + 'e' + (value[1] ? (+value[1] - exp) : -exp));
}

// JS equivalent of PHP's $_GET
function parseURLParams(url) {
    var queryStart = url.indexOf("?") + 1,
        queryEnd   = url.indexOf("#") + 1 || url.length + 1,
        query = url.slice(queryStart, queryEnd - 1),
        pairs = query.replace(/\+/g, " ").split("&"),
        parms = {}, i, n, v, nv;

    if (query === url || query === "") {
        return;
    }

    for (i = 0; i < pairs.length; i++) {
        nv = pairs[i].split("=");
        n = decodeURIComponent(nv[0]);
        v = decodeURIComponent(nv[1]);

        if (!parms.hasOwnProperty(n)) {
            parms[n] = [];
        }

        parms[n].push(nv.length === 2 ? v : null);
    }
    return parms;
}

function flattenJSON(object) {
  var x = new Array();
  for (var k in object) {
    x.push(object[k]);
  }

  return x;
}