listener["initialize"] = function()
{
  listener.available_commands = ["startserver","help","ls","clear"];
  listener.help_commands = ["start the websocket server", "what this did","unsupported...","clear the terminal"];
}

isfile = function(text)
{
  if (text==".." || text=="../") return false;
  text = text.split(".");
  if (text.length>1) return true;
  return false;
}

listener["parse_command"] = function(msg)
{
  var out = {text: [], color: ""};

  data = msg;
  data = data.split(" ");
  msg = data[0];

  // create a new line like a terminal when enter is hit
  if (msg=="") return out;

  if (msg=="send")
  {
  }
  else if (msg=="startserver")
  {
    startserver();
  }
  else if (msg=="killserver")
  {
    killserver();
  }
  else if (msg=="help")
  {
    out = print_help(out);
  }
  else if (msg=="ls")
  {
    wv.socketUt.send("ls|");

    // the server will send a response back to wv and we can get the values of the directory
    var terminal = document.getElementById("brframe");
    terminal.listener.print_directory();
  }
  else if (msg=="cd")
  {
    if (data.length==1)
    {}
    else if (isfile(data[1]))
      alert("cannot cd into this!");
    else
      wv.socketUt.send("cd|"+data[1]);
  }
  else if (msg=="plot")
  {
    // need to do a lot of work to parse the rest of the command
    parse_plot(data);
  }
  else
  {
  }

  return out;
}

parse_plot = function(command)
{
  print_message(command,"yellow",false,false);
  wv.socketUt.send("plot|"+command[1]); // for now...
}

listener["print_directory"] = function()
{
  text = "implement!";
  print_message( text, "white" , false , false );
}


startserver = function()
{
  print_message( "implement!" );
}

killserver = function()
{
  print_message( "implement!" );
}

print_help = function(out)
{
  for (i=0;i<listener.available_commands.length;i++)
  {
    out.text[i] = listener.available_commands[i]+": "+listener.help_commands[i];
  }
  out.color = "yellow";
  return out;
}
