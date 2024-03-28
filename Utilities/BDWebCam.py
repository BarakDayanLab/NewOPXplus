import os.path

import cv2
import time


class BDWebCam:

    def __init__(self, folder, device_index=0):
        self.device_index = device_index  # the camera device index (they may be more than one connected
        self.webcam = None
        self.folder = folder
        pass

    def connect(self):
        # Try either VideoCapture(0) or VideoCapture(1) based on your camera availability
        self.webcam = cv2.VideoCapture(self.device_index)
        pass

    def disconnect(self):
        if self.webcam is not None:
            self.webcam.release()
        self.webcam = None
        pass

    def start_live_stream(self, save_images, delay):

        last_snapshot = int(time.time())

        while True:
            # Capture the video frame by frame
            ret, frame = self.webcam.read()

            # Display the resulting frame
            cv2.imshow('Live Stream (press q to quit)', frame)

            # Check if we should save the
            if save_images:
                now_epoch = int(time.time())
                if (now_epoch - last_snapshot) > delay:
                    last_snapshot = now_epoch
                    self.save_image(frame)

            # the 'q' button is set as the
            # quitting button you may use any
            # desired button of your choice
            if cv2.waitKey(1) & 0xFF == ord('q'):
                break
        pass

    def snap(self):
        check, frame = self.webcam.read()
        return check, frame

    def timelapse(self, delay=1, number_of_frames=-1):
        """
        Takes consequent images from webcam.

        delay - delay between frames (in N seconds)
        number_of_frames - how many images to take (-1 for endless)
        folder - where to save the images

        """
        count = 0
        while count < number_of_frames:
            # Capture the video frame by frame
            ret, frame = self.webcam.read()

            # Save image to folder
            self.save_image(image=frame)

            count += 1
            time.sleep(delay)

        pass

    def save_image(self, image):
        # Save image to folder
        current_date_time = time.strftime("%Y%m%d_%H%M%S")
        file_name = f'{current_date_time}.jpg'
        full_path = os.path.join(self.folder, file_name)
        cv2.imwrite(filename=full_path, img=image)
        print(f'Took snapshot - {full_path}')

    @staticmethod
    def test():

        save_path = r'C:\temp\GaussMeter'

        webcam = BDWebCam(folder=save_path)
        webcam.connect()

        mode = 'LIVE'

        if mode == 'LIVE':
            webcam.start_live_stream(save_images=True, delay=60*5)
        elif mode == 'TIMELAPSE':
            webcam.timelapse(delay=60*5, number_of_frames=-1)  # Every 5 minutes
        elif mode == 'SNAP':
            check, frame = webcam.snap()
            cv2.imshow("Image", frame)
            cv2.waitKey(0)

        webcam.disconnect()

        pass


        #cv2.imwrite(filename=r'<Your Directory>\saved_img.jpg', img=frame)

        webcam.release()
        pass

if __name__ == "__main__":

    BDWebCam.test()
    pass