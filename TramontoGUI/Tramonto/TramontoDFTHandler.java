package Tramonto;
import RemoteFile.*;
import XML.*;
import IdeaComm.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.*;
import java.io.*;

public class TramontoDFTHandler 
    extends MessageHandler 
    implements ItemListener, ActionListener
{

    // Swing components used for interacting with the user
    JDialog theDialogWindow;
    Container theDialogContent;
    JButton theOkButton;
    JButton theCancelButton;
    JCheckBox theValidateCheckBox;
    JCheckBox theGenInputFileCheckBox;
    JCheckBox theValidateSurfacesFileCheckBox;
    /*
    ButtonGroup theSimStartGroup;
    JRadioButton theDoNotStartSimRadioButton;
    JRadioButton theStartLocalSimRadioButton;
    JRadioButton theStartRemoteSimRadioButton;
    */

    // XML object passed to the handler when it is created
    XMLObject theInputParameterMessage;

    // parser used to manipulate the input message
    TramontoXMLParser theParser;

    //
    // Default constructor
    //
    public TramontoDFTHandler() {
	// set up dialog box to be shown later
	setupProcessingDialog();
    }


    //
    // entry point for calls from Maui/Idea main program
    //
    public void handleTramontoDFT(IdeaMessage msg) {
	// save the message
	theInputParameterMessage = msg.getBody();

	// pop up the dialog window
        theDialogWindow.setVisible(true);
    }


    //
    // ActionListerner interface routine to catch button events
    //
    public void actionPerformed(ActionEvent evt) {
	// do the work if theOkButton has been clicked
	if( evt.getSource() == theOkButton ) {
	    if(theValidateCheckBox.isSelected() == true) {
		boolean validInput = validateInput();
	        if( validInput==true ) {
		    if( theGenInputFileCheckBox.isSelected() == true ) {
			generateInputFile();
		    }
		    if( theValidateSurfacesFileCheckBox.isSelected() == true ) {
			checkSurfacesFile();
		    }

		    /*
		    if((makeFileResult == true) &&
		       (theStartLocalSimRadioButton.isSelected() == true)) {
			startLocalSimulation();
		    }
		    else if(theStartRemoteSimRadioButton.isSelected() == true) {
			startRemoteSimulation();
		    }
		    */
		}
	    }
	    // did the work so now close the window.
	    theDialogWindow.dispose();
	}
	// if theCancelButton was pressed then destroy the window
	else if( evt.getSource() == theCancelButton ) {
	    theDialogWindow.dispose();
	}
    }


    //
    // ItemListener interface routine to catch checkbox events
    //
    public void itemStateChanged(ItemEvent evt) {
	if( evt.getSource() == theValidateCheckBox ) {
	    if(theValidateCheckBox.isSelected()==true) {
		theGenInputFileCheckBox.setEnabled(true);
		theValidateSurfacesFileCheckBox.setEnabled(true);
		/*
		if(theGenInputFileCheckBox.isSelected()==true) {
		    theDoNotStartSimRadioButton.setEnabled(true);
		    theStartLocalSimRadioButton.setEnabled(true);
		    // Remote simulations are disabled until the 
		    // IdeaServer/Client is debuged
		    theStartRemoteSimRadioButton.setEnabled(false);
		}
		*/
	    }
	    else {
		theGenInputFileCheckBox.setEnabled(false);
		theValidateSurfacesFileCheckBox.setEnabled(false);
		/*
                theDoNotStartSimRadioButton.setEnabled(false);
		theStartLocalSimRadioButton.setEnabled(false);
		theStartRemoteSimRadioButton.setEnabled(false);
		*/
	    }
	}
	else if( evt.getSource() == theGenInputFileCheckBox ) {
	    /*
	    if( theGenInputFileCheckBox.isSelected() == true) {
                theDoNotStartSimRadioButton.setEnabled(true);
		theStartLocalSimRadioButton.setEnabled(true);
		// Remote simulations are disabled until the 
		// IdeaServer/Client is debuged
		theStartRemoteSimRadioButton.setEnabled(false);
	    }
	    else {
                theDoNotStartSimRadioButton.setEnabled(false);
		theStartLocalSimRadioButton.setEnabled(false);
		theStartRemoteSimRadioButton.setEnabled(false);
	    }
	    */
	}

		
    }

    //
    // Private routine to configure the dialog box used to get
    // the user's direction for the type of processing s/he 
    // desires. 
    //
    private void setupProcessingDialog() {
        // set up dialog window to find out what the user wants done
        theDialogWindow = new JDialog();
        theDialogWindow.setTitle("TramontoDFT Processing");
        theDialogWindow.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        theDialogContent = theDialogWindow.getContentPane();

	// make space for buttons on the bottom of the dialogbox
        JPanel theButtonPanel = new JPanel();

	// make the standard "ok" and "cancel" buttons, 
	// let this object listen to the buttons
	// and add them to the panel
        theOkButton = new JButton("Ok");
	theOkButton.addActionListener(this);
        theButtonPanel.add(theOkButton);
        theCancelButton = new JButton("Cancel");
	theCancelButton.addActionListener(this);
        theButtonPanel.add(theCancelButton);
        theDialogContent.add(theButtonPanel, BorderLayout.SOUTH);

        // processing options for the XML message sent to this object
	// Validating input parameters is a prerequisite for 
	// creating a simulation input file.  
	// creating a simulation input file is a prerequisite for
	// starting a simulation on the local or a remote machine.

	// a panel to hold our checkboxes and radiobuttons
        JPanel theSelectionPanel = new JPanel();

	// set the layout to infinite rows, one column
        theSelectionPanel.setLayout(new GridLayout(0,1));
        theSelectionPanel.add( new JLabel("Please select desired processing options"));
	theValidateCheckBox = new JCheckBox("Validate input parameters");
        theValidateCheckBox.setSelected(true);
	// set up this object to listen to the checkbox's This is done so that we can
	// enable and disable the appropriate options depending on what's checked
	theValidateCheckBox.addItemListener(this);
        theGenInputFileCheckBox = new JCheckBox("Generate simulation input file");
	theGenInputFileCheckBox.setSelected(true);
        theValidateSurfacesFileCheckBox = new JCheckBox("Check dft_surfaces.dat file for errors");
        theValidateSurfacesFileCheckBox.setSelected(true);

	// set up this object to listen to the checkbox's 
	/*
        theGenInputFileCheckBox.addItemListener(this);
        theSimStartGroup = new ButtonGroup();
        theDoNotStartSimRadioButton = new JRadioButton("Do not start simulation");
	theDoNotStartSimRadioButton.setSelected(true);
        theSimStartGroup.add( theDoNotStartSimRadioButton);
        theStartLocalSimRadioButton = new JRadioButton("Start simulation on local machine");
        theSimStartGroup.add( theStartLocalSimRadioButton );
        theStartRemoteSimRadioButton =
	    new JRadioButton("Start simulation on remote machine");
	*/
	// Remote simulations are disabled until the IdeaServer/Client is debuged
	/*
	theStartRemoteSimRadioButton.setEnabled(false);
	theSimStartGroup.add( theStartRemoteSimRadioButton );
	*/

        // build up the panel
        theSelectionPanel.add( theValidateCheckBox );
        theSelectionPanel.add( theGenInputFileCheckBox );
        theSelectionPanel.add( theValidateSurfacesFileCheckBox );
        //theSelectionPanel.add( theDoNotStartSimRadioButton );
        //theSelectionPanel.add( theStartLocalSimRadioButton );
        //theSelectionPanel.add( theStartRemoteSimRadioButton );
        theDialogContent.add( theSelectionPanel, BorderLayout.CENTER );

	// pack up the window and place it so it can be sized and centered 
        theDialogWindow.pack();

	// center the dialogbox
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
	Dimension windowSize = theDialogWindow.getSize();
        int xPos = (int) ((screenSize.getWidth() - windowSize.getWidth()) / 2);
        int yPos = (int) ((screenSize.getHeight() - windowSize.getHeight()) / 2);
        theDialogWindow.setLocation(xPos, yPos);

    }


    //
    // Try and validate the input contained in the XML message passed to 
    // this object from Maui/Idea.  Return true if XML passed as valid
    // otherwise return false.
    //
    private boolean validateInput() {
	theParser = new TramontoXMLParser(theInputParameterMessage);
	try {
	    theParser.validateXML();
	}
	catch(RuntimeException rte ) {
	    // get the validation error message
	    String theMessage = rte.getMessage();

	    String theLocation;
	    StringWriter theError = new StringWriter();
	    PrintWriter thePrintWriter = new PrintWriter( theError );
	    // get the stacktrace
	    rte.printStackTrace( thePrintWriter );
	    // extract the location
	    String theFullTrace = theError.toString();
	    int startLoc = theFullTrace.indexOf("\n\tat ");
	    // look over package and class names
	    startLoc = theFullTrace.indexOf(".", startLoc) + 1;
	    startLoc = theFullTrace.indexOf(".", startLoc) + 1;
	    // chop off routine arguments
	    int endLoc = theFullTrace.indexOf("(", startLoc);
	    theLocation = theFullTrace.substring(startLoc, endLoc);
	    theMessage = theMessage + "\nIn routine: " + theLocation;

	    // Throw up an error box
	    JOptionPane.showMessageDialog(null, theMessage, "Input Parameter Validation Error", JOptionPane.ERROR_MESSAGE);
	    return false;
	}
	return true;
    }


    //
    // Generate an input file for the simulation.
    // Return true if file was created, otherwise false
    //
    private boolean generateInputFile() {
	return theParser.generateInputFile();
    }

    //
    // check the dft_surfaces file
    //
    private boolean checkSurfacesFile() {
	try {
	    theParser.validateSurfacesFile();
	}
	catch(RuntimeException rte ) {
	    // get the validation error message
	    String theMessage = rte.getMessage();

	    // Throw up an error box
	    JOptionPane.showMessageDialog(null, theMessage, "Surfaces file Error:", JOptionPane.ERROR_MESSAGE);
	    return false;
	}
	return true;
    }


    //
    // Start simulation on local host
    //
    private void startLocalSimulation() {
	System.out.println("\tIn startLocalSimulation()");
    }


    //
    // Start simulation on remote host
    //
    private void startRemoteSimulation() {
	System.out.println("\tIn startRemoteSimulation()");
    }
}
