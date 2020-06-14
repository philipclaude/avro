/* utilities */
if( !String.prototype.trim ) {
  String.prototype.trim = function () {
    return this.replace(/^\s+|\s+$/g,'');
  };
}

/* parsing callback */
var parse_command = function(text) {
	text = text.trim();
  term = document.getElementById("brframe");

  // catch a "clear"
  if (text=="clear")
  {
    term.clear();
    return;
  }

  out = term.listener.parse_command(text);
  if (out.text!="")
  {
    for (i=0;i<out.text.length;i++)
      term.print_message(out.text[i],out.color,false);
  }
}

/* initialize a new terminal */
var InitTerminal = function() {
  t = new Terminal();
  document.body.appendChild(t.html)
  create_line();
  return t;
}

/* create new terminal line */
var create_line = function() {

  var term = document.getElementById("brframe");

  /* the new input line */
  input_ = document.createElement('p');

  /* children of the input line */
  prefix_  = document.createElement('span');
  command_ = document.createElement('input');

  command_.style.background = 'black';
  command_.style.color = 'white';
  command_.style.border = 'none';
  command_.style.outline = 'none';
  command_.style.width = "95%";//term.clientWidth -30;
  command_.style.zindex = -100;
  command_.style.fontFamily = "monospace";
  command_.style.fontWeight = "bold";
  command_.style.spellCheck = "false";
  command_.style.maxlength = 250;

  //command_.style.fontSize = "10px";

  prefix_.innerHTML = '<b> >> </b>';
  prefix_.style.background = 'black';
  prefix_.style.color = 'white';
  prefix_.style.fontWeight = "bold";
  prefix_.style.fontFamily = "monospace";
  //prefix_.style.fontSize = "10px";

  input_.appendChild(prefix_);
  input_.appendChild(command_);

  input_.style.width = term.clientWidth -30;

  if (!term.firstTime)
  {
	   command_.focus();
  }
  else
  {
	   term.firstTime = false;
  }

  //input_.style.height = "16px";
  //input_.style.lineHeight = "14px";
  input_.style.margin = '0px 0px 0px 0px';

  term.appendChild(input_);

  input_.scrollIntoView();

  term.activeCommand = command_;

  term.activeCommand.onfocus = function() { wv.terminalON = true; };
  term.activeCommand.onblur = function() { wv.terminalON = false; };

  // key up function
  command_.onkeyup = function(e)
  {

    if (e.which===13)
    {
      // enter was pressed

      term.commandlist.push(command_);

  	  e.preventDefault();
  	  command_.readOnly = true;
  	  parse_command(command_.value);
  	  create_line();
  	  setTimeout(command_.focus(),0);
  	}
  }

  // post message function
  term.print_message = function( msg , col , newline )
  {

    msg = msg.split('\n');
    for (var i=0;i<msg.length;i++)
    {

      /* the new input line */
    	input_ = document.createElement('p');

    	/* children of the input line */
      msg_ = document.createElement('span');

    	msg_.innerHTML = '<b>  </b>' +msg[i];
    	msg_.style.background = 'black';
    	msg_.style.color = col;
      msg_.style.fontWeight = "bold";
      msg_.style.fontFamily = "monospace";
      //msg_.style.fontSize = "10px";
      //msg_.style.lineHeight = "14px";

      //input_.style.height = "10px";
      //input_.style.lineHeight = "10px";
      input_.style.margin = '2px 0px 2px 0px';

      input_.appendChild(msg_);
      this.appendChild(input_);
    }

    if (newline)
      create_line();
  }

  term.set_listener = function( _listener )
  {
    this.listener = _listener;
  }

  term.clear = function()
  {
    while (this.childNodes.length!=0)
    {
      children = this.childNodes;
      this.removeChild(children[0]);
    }
  }
}

var Terminal = ( function () {

	var makeTerminal = function () {

		this.html = document.createElement('brframe')
		this.html.id = "brframe";

		this.html.firstTime = true;

		this.html.style.background = 'black';
		this.html.style.color = 'white';
		//this.html.style.fontSize = '8pt';
		this.html.style.width = '100%';
		this.html.style.border = '2px solid #ccc';
		this.html.style.padding = '2px';
		this.html.style.float = 'left';
		this.html.style.margin = '10px 0px 10px 0px';
		this.html.style.height = '200px';
		this.html.style.overflowY = "scroll";
		this.html.tabIndex = 0;

		this.html.style.fontFamily = 'Monaco, Courier';

		this.html.commandlist = [];
		this.html.outputlist  = [];

		this.html.onfocus = function() {
      this.activeCommand.focus();
      //wv.terminalON = true;

		}

    this.html.onblur = function() {
      //wv.terminalON = false;
    }

	}

	return makeTerminal;

}())
