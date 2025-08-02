import classes from "./Contact.module.css";
import ContactForm from "./ContactForm";
import Dialog from "./Dialog";

const ContactLink = () => (
  <Dialog
    title="Contact Us"
    trigger={
      <button className={classes.link} type="button" data-tooltip="Contact us">
        Contact
      </button>
    }
  >
    <ContactForm />
  </Dialog>
);

export default ContactLink;
