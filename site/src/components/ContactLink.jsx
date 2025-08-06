import Button from "./Button";
import classes from "./Contact.module.css";
import ContactForm from "./ContactForm";
import Dialog from "./Dialog";

const ContactLink = () => (
  <Dialog
    title="Contact Us"
    trigger={
      <Button
        id="footer-contact"
        className={classes.link}
        data-tooltip="Contact us"
      >
        Contact
      </Button>
    }
  >
    <ContactForm />
  </Dialog>
);

export default ContactLink;
