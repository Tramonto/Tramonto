package Tramonto;
import RemoteFile.*;
import XML.*;
import IdeaComm.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.*;
import java.io.*;

public class TramontoStateParametersHandler 
    extends MessageHandler 
{

    // XML object passed to the handler when it is created
    XMLObject theInputParameterMessage;

    // parser used to manipulate the input message
    TramontoXMLParser theParser;

    //
    // Default constructor
    //
    public TramontoStateParametersHandler() {
    }


    //
    // entry point for calls from Maui/Idea main program
    //
    public void handleTramontoStateParameters(IdeaMessage msg) {
	// save the message
	theInputParameterMessage = msg.toXMLObject();
	validateInput();
    }



    //
    // Try and validate the input contained in the XML message passed to 
    // this object from Maui/Idea.  Return true if XML passed is valid
    // otherwise return false.
    //
    private boolean validateInput() {
	theParser = new TramontoXMLParser(theInputParameterMessage);
	try {
	    theParser.validateFunctionalSwitchParameters();
	    theParser.validatePotentialTypeParameters();
	    theParser.validateStatePointParameters();
	    
	}
	catch(RuntimeException rte ) {
	    // get the validation error message
	    String theMessage = rte.getMessage();

	    // Throw up an error box
	    JOptionPane.showMessageDialog(null, theMessage, "Input Parameter Validation Error", JOptionPane.ERROR_MESSAGE);
	    return false;
	}
	return true;
    }
}
